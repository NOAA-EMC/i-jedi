#include <netcdf.h>

#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/mom6/IoMOM6.h"
#include "ijedi/Geometry/mom6/GeometryMOM6.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoMOM6> makerIoMOM6_("mom6");
  // -------------------------------------------------------------------------------------------------
  IoMOM6::IoMOM6(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()),
        geom_(geom),
        datapath_(params.datapath),
        filename_(params.filename.value().value_or(""))
  {
    util::Timer timer(classname(), "IoMOM6");
    oops::Log::trace() << classname() << " constructor starting" << std::endl;
    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  IoMOM6::~IoMOM6()
  {
    util::Timer timer(classname(), "~IoMOM6");
    oops::Log::trace() << classname() << " destructor starting" << std::endl;
    oops::Log::trace() << classname() << " destructor done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoMOM6::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                    const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read state");
    oops::Log::trace() << classname() << " read state starting" << std::endl;

    // Build file path: if filename_ is absolute use it directly, else prepend datapath_
    const std::string filepath = (!filename_.empty() && filename_[0] == '/')
        ? filename_
        : datapath_ + "/" + filename_;

    // Open netCDF file (read-only, all ranks)
    int ncid;
    ASSERT_MSG(nc_open(filepath.c_str(), NC_NOWRITE, &ncid) == NC_NOERR,
               classname() + "::read: cannot open " + filepath);

    // Read horizontal and vertical file dimensions.
    // Convention (same as GeometryMOM6): global_index = jG * lonh + iG + 1
    int dimid;
    size_t file_ni, file_nj, file_nz;
    ASSERT_MSG(nc_inq_dimid(ncid, "lonh",  &dimid) == NC_NOERR,
               classname() + "::read: missing dimension 'lonh' in " + filepath);
    nc_inq_dimlen(ncid, dimid, &file_ni);
    ASSERT_MSG(nc_inq_dimid(ncid, "lath",  &dimid) == NC_NOERR,
               classname() + "::read: missing dimension 'lath' in " + filepath);
    nc_inq_dimlen(ncid, dimid, &file_nj);
    ASSERT_MSG(nc_inq_dimid(ncid, "Layer", &dimid) == NC_NOERR,
               classname() + "::read: missing dimension 'Layer' in " + filepath);
    nc_inq_dimlen(ncid, dimid, &file_nz);

    // Sanity-check vertical dimension against geometry
    const auto & geomMOM6read = dynamic_cast<const GeometryMOM6 &>(geom_.geometryImpl());
    ASSERT_MSG(static_cast<int>(file_nz) == geomMOM6read.numLevels(),
               classname() + "::read: file has " + std::to_string(file_nz)
               + " layers but geometry has " + std::to_string(geomMOM6read.numLevels()));

    // Default MOM6 variable name mapping: JEDI name → file variable name.
    // U/V are on staggered grids and are not handled here yet.
    static const std::unordered_map<std::string, std::string> kDefaultNames = {
      {"sea_water_potential_temperature",   "Temp"},
      {"sea_water_conservative_temperature","Temp"},
      {"sea_water_practical_salinity",      "Salt"},
      {"sea_water_salinity",                "Salt"},
      {"sea_water_cell_thickness",          "h"},
      {"sea_surface_height_above_geoid",    "ave_ssh"},
    };

    // Global index view: g = jG * file_ni + iG + 1  (same scheme in both
    // mom6FunctionSpace_ and functionSpace_; see GeometryMOM6)
    const auto & fs =
        atlas::functionspace::NodeColumns(x.field(0).functionspace());
    auto g_view = atlas::array::make_view<atlas::gidx_t, 1>(
        fs.mesh().nodes().global_index());
    const int nNodes = static_cast<int>(fs.mesh().nodes().size());

    for (auto & field : x) {
      const std::string jediName = field.name();
      const int nLevels = field.shape(1);
      const bool is3D   = (nLevels > 1);

      // Resolve file variable name
      std::string fileVarName;
      if (fileionames.has(jediName)) {
        fileVarName = fileionames.getString(jediName);
      } else {
        const auto it = kDefaultNames.find(jediName);
        if (it == kDefaultNames.end()) {
          oops::Log::warning() << classname() << "::read: no MOM6 file mapping for '"
                               << jediName << "', skipping" << std::endl;
          continue;
        }
        fileVarName = it->second;
      }

      int varid;
      if (nc_inq_varid(ncid, fileVarName.c_str(), &varid) != NC_NOERR) {
        oops::Log::warning() << classname() << "::read: variable '" << fileVarName
                             << "' not found in file, skipping '" << jediName
                             << "'" << std::endl;
        continue;
      }

      auto view = atlas::array::make_view<double, 2>(field);

      if (is3D) {
        // Read full (1, nz, nj, ni) slab at Time=0
        std::vector<double> buf(file_nz * file_nj * file_ni);
        const std::vector<size_t> start = {0, 0, 0, 0};
        const std::vector<size_t> count = {1, file_nz, file_nj, file_ni};
        nc_get_vara_double(ncid, varid, start.data(), count.data(), buf.data());

        if (fileioscaling.has(jediName)) {
          const double scale = fileioscaling.getDouble(jediName);
          for (auto & v : buf) v *= scale;
        }

        for (int n = 0; n < nNodes; ++n) {
          const int iG = static_cast<int>((g_view(n) - 1) % static_cast<int>(file_ni));
          const int jG = static_cast<int>((g_view(n) - 1) / static_cast<int>(file_ni));
          for (int k = 0; k < nLevels; ++k) {
            view(n, k) = buf[static_cast<size_t>(k)  * file_nj * file_ni
                           + static_cast<size_t>(jG) * file_ni
                           + static_cast<size_t>(iG)];
          }
        }

      } else {
        // 2D variable: read full (1, nj, ni) slab at Time=0
        std::vector<double> buf(file_nj * file_ni);
        const std::vector<size_t> start = {0, 0, 0};
        const std::vector<size_t> count = {1, file_nj, file_ni};
        nc_get_vara_double(ncid, varid, start.data(), count.data(), buf.data());

        if (fileioscaling.has(jediName)) {
          const double scale = fileioscaling.getDouble(jediName);
          for (auto & v : buf) v *= scale;
        }

        for (int n = 0; n < nNodes; ++n) {
          const int iG = static_cast<int>((g_view(n) - 1) % static_cast<int>(file_ni));
          const int jG = static_cast<int>((g_view(n) - 1) / static_cast<int>(file_ni));
          view(n, 0) = buf[static_cast<size_t>(jG) * file_ni + static_cast<size_t>(iG)];
        }
      }

      oops::Log::info() << classname() << "::read: " << jediName
                        << " <- " << fileVarName << " from " << filepath << std::endl;
    }

    nc_close(ncid);

    // Write to netcdf using the FieldSet helpers for debugging purposes
    // TODO: remove this once read is verified and tested, or make it conditional on a debug flag
    oops::Log::trace() << "write MOM6 (fieldset) starting" << std::endl;
    {
      std::string debugPath = filepath;
      const std::string ext = ".nc";
      if (debugPath.size() > ext.size() &&
          debugPath.compare(debugPath.size() - ext.size(), ext.size(), ext) == 0)
        debugPath.erase(debugPath.size() - ext.size());
      debugPath += "_debug";
      eckit::LocalConfiguration debugConf;
      debugConf.set("filepath", debugPath);
      util::writeFieldSet(geom_.getComm(), debugConf, x);
    }
    oops::Log::trace() << "write MOM6 (fieldset) done" << std::endl;

    // Halo exchange to fill ghost nodes from their owning ranks
    for (auto & field : x) {
      atlas::functionspace::NodeColumns(field.functionspace()).haloExchange(field);
    }

    oops::Log::trace() << classname() << " read state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoMOM6::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                     const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;

    // Build file path: if filename_ is absolute use it directly, else prepend datapath_
    const std::string filepath = (!filename_.empty() && filename_[0] == '/')
        ? filename_
        : datapath_ + "/" + filename_;

    const auto & comm = geom_.getComm();
    const bool isRoot = (comm.rank() == 0);

    // Get geometry info
    const auto & geomMOM6 = dynamic_cast<const GeometryMOM6 &>(geom_.geometryImpl());
    const int ni = geomMOM6.niEff();
    const int nj = geomMOM6.njEff();
    const int nz = geomMOM6.numLevels();

    // Default MOM6 variable name mapping
    static const std::unordered_map<std::string, std::string> kDefaultNames = {
      {"sea_water_potential_temperature",   "Temp"},
      {"sea_water_conservative_temperature","Temp"},
      {"sea_water_practical_salinity",      "Salt"},
      {"sea_water_salinity",                "Salt"},
      {"sea_water_cell_thickness",          "h"},
      {"sea_surface_height_above_geoid",    "ave_ssh"},
    };

    // Resolve variable names for all fields (needed on all ranks)
    std::vector<std::string> fileVarNames;
    for (const auto & field : x) {
      const std::string jediName = field.name();
      std::string fileVarName;
      if (fileionames.has(jediName)) {
        fileVarName = fileionames.getString(jediName);
      } else {
        const auto it = kDefaultNames.find(jediName);
        if (it == kDefaultNames.end()) { fileVarName = ""; }
        else { fileVarName = it->second; }
      }
      fileVarNames.push_back(fileVarName);
    }

    // Rank 0 creates the file and defines all dimensions and variables
    int ncid = -1;
    std::unordered_map<std::string, int> varids;
    if (isRoot) {
      std::remove(filepath.c_str());
      ASSERT_MSG(nc_create(filepath.c_str(), NC_CLOBBER, &ncid) == NC_NOERR,
                 classname() + "::write: cannot create " + filepath);

      int dim_time, dim_layer, dim_lath, dim_lonh;
      nc_def_dim(ncid, "Time",  NC_UNLIMITED, &dim_time);
      nc_def_dim(ncid, "Layer", nz,           &dim_layer);
      nc_def_dim(ncid, "lath",  nj,           &dim_lath);
      nc_def_dim(ncid, "lonh",  ni,           &dim_lonh);

      int fi = 0;
      for (const auto & field : x) {
        const std::string & fileVarName = fileVarNames[fi++];
        if (fileVarName.empty()) continue;
        if (varids.count(fileVarName)) continue;  // already defined (e.g. two JEDI names → same var)
        int varid;
        if (field.shape(1) > 1) {
          int dimids[4] = {dim_time, dim_layer, dim_lath, dim_lonh};
          nc_def_var(ncid, fileVarName.c_str(), NC_DOUBLE, 4, dimids, &varid);
        } else {
          int dimids[3] = {dim_time, dim_lath, dim_lonh};
          nc_def_var(ncid, fileVarName.c_str(), NC_DOUBLE, 3, dimids, &varid);
        }
        varids[fileVarName] = varid;
      }
      nc_enddef(ncid);
    }

    // Global index and ghost views (all ranks)
    const auto & fs =
        atlas::functionspace::NodeColumns(x.field(0).functionspace());
    auto g_view = atlas::array::make_view<atlas::gidx_t, 1>(
        fs.mesh().nodes().global_index());
    auto ghost_view = atlas::array::make_view<int, 1>(
        fs.mesh().nodes().ghost());
    const int nNodes = static_cast<int>(fs.mesh().nodes().size());

    // Write each field: each rank fills only its owned (non-ghost) nodes into
    // a global buffer, then reduce(sum) to root which writes to file.
    int fi = 0;
    for (const auto & field : x) {
      const std::string jediName  = field.name();
      const std::string & fileVarName = fileVarNames[fi++];
      if (fileVarName.empty()) continue;

      auto view = atlas::array::make_view<double, 2>(field);
      const int nLevels = field.shape(1);

      if (nLevels > 1) {
        std::vector<double> buf(static_cast<size_t>(nz * nj * ni), 0.0);
        for (int n = 0; n < nNodes; ++n) {
          if (ghost_view(n)) continue;   // skip ghost nodes
          const int iG = static_cast<int>((g_view(n) - 1) % ni);
          const int jG = static_cast<int>((g_view(n) - 1) / ni);
          for (int k = 0; k < nLevels; ++k) {
            const double val = view(n, k);
            if (std::abs(val) < 1.0e30)
              buf[static_cast<size_t>(k * nj * ni + jG * ni + iG)] = val;
          }
        }
        comm.reduceInPlace(buf.data(), buf.size(), eckit::mpi::sum(), 0);
        if (isRoot) {
          const size_t start[4] = {0, 0, 0, 0};
          const size_t count[4] = {1, static_cast<size_t>(nz),
                                      static_cast<size_t>(nj),
                                      static_cast<size_t>(ni)};
          nc_put_vara_double(ncid, varids.at(fileVarName), start, count, buf.data());
        }
      } else {
        std::vector<double> buf(static_cast<size_t>(nj * ni), 0.0);
        for (int n = 0; n < nNodes; ++n) {
          if (ghost_view(n)) continue;   // skip ghost nodes
          const int iG = static_cast<int>((g_view(n) - 1) % ni);
          const int jG = static_cast<int>((g_view(n) - 1) / ni);
          const double val = view(n, 0);
          if (std::abs(val) < 1.0e30)
            buf[static_cast<size_t>(jG * ni + iG)] = val;
        }
        comm.reduceInPlace(buf.data(), buf.size(), eckit::mpi::sum(), 0);
        if (isRoot) {
          const size_t start[3] = {0, 0, 0};
          const size_t count[3] = {1, static_cast<size_t>(nj), static_cast<size_t>(ni)};
          nc_put_vara_double(ncid, varids.at(fileVarName), start, count, buf.data());
        }
      }
      oops::Log::info() << classname() << "::write: " << jediName
                        << " -> " << fileVarName << " to " << filepath << std::endl;
    }

    if (isRoot) nc_close(ncid);
    oops::Log::trace() << classname() << " write state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoMOM6::print(std::ostream &os) const
  {
    os << classname() << " Io for MOM6 restarts and histories";
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

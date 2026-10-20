#include "ijedi/Io/mom6/WriteMOM6Netcdf.h"

#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/mom6/GeometryMOM6.h"

namespace ijedi {

namespace {
void ncCheck(int status, const std::string & what, const std::string & filepath) {
  if (status != NC_NOERR)
    throw eckit::Exception("writeMOM6Netcdf: " + what + " failed for " + filepath + ": "
                           + nc_strerror(status), Here());
}
}  // namespace

// ---------------------------------------------------------------------------
void writeMOM6Netcdf(const std::string & filepath,
                     const atlas::FieldSet & x,
                     const std::vector<std::string> & fileVarNames,
                     const int ni, const int nj, const int nz,
                     const eckit::mpi::Comm & comm,
                     const bool nativeVelocity,
                     const mom6::VelocityGridInfo * velGrid) {
  const bool isRoot = (comm.rank() == 0);

  // Velocity pairs to write on the native MOM6 C-grid (grid-relative, u on
  // east faces, v on north faces).  eastIdx -> northIdx.
  std::unordered_map<int, int> nativePairs;
  std::unordered_set<int> nativeNorth;
  if (nativeVelocity) {
    for (const auto & pair : mom6::velocityPairs()) {
      int iE = -1, iN = -1;
      for (int f = 0; f < x.size(); ++f) {
        if (fileVarNames[f].empty()) continue;
        if (x.field(f).name() == pair.first)  iE = f;
        if (x.field(f).name() == pair.second) iN = f;
      }
      if (iE < 0 && iN < 0) continue;
      if (iE < 0 || iN < 0)
        throw eckit::UserError(
            "writeMOM6Netcdf: velocity grid native needs both " + pair.first +
            " and " + pair.second, Here());
      nativePairs[iE] = iN;
      nativeNorth.insert(iN);
    }
    ASSERT_MSG(nativePairs.empty() || !isRoot || velGrid != nullptr,
               "writeMOM6Netcdf: native velocity requires grid rotation info");
  }

  // --- 1. Gather every field to rank 0 using Atlas NodeColumns::gather ---
  const auto fs =
      atlas::functionspace::NodeColumns(x.field(0).functionspace());

  // The JEDI function space is unstructured and only contains active points
  // (ocean + fringe); land points are omitted.  After fs.gather() the global
  // field has nActiveGlobal entries -- NOT ni*nj.  We therefore also gather
  // the global-index field so we know where each gathered node sits on the
  // full structured grid.
  //
  // Convention: global_index = jG * ni + iG + 1  (1-based)

  // 1a. Build a local field holding each node's global index as double,
  //     then gather it to root.
  auto localGidx = atlas::array::make_view<atlas::gidx_t, 1>(
      fs.nodes().global_index());
  atlas::Field gidxField = fs.createField<double>(
      atlas::option::name("__gidx") | atlas::option::levels(1));
  {
    auto gv = atlas::array::make_view<double, 2>(gidxField);
    for (atlas::idx_t n = 0; n < fs.nodes().size(); ++n)
      gv(n, 0) = static_cast<double>(localGidx(n));
  }
  atlas::Field globalGidxField = fs.createField<double>(
      atlas::option::name("__gidx") | atlas::option::levels(1)
      | atlas::option::global());
  fs.gather(gidxField, globalGidxField);

  // 1b. Gather all data fields.
  std::vector<atlas::Field> globalFields;
  globalFields.reserve(x.size());

  int fi = 0;
  for (const auto & field : x) {
    const std::string & varName = fileVarNames[fi++];
    if (varName.empty()) {
      globalFields.emplace_back();          // placeholder
      continue;
    }
    const int nLevels = field.shape(1);
    atlas::Field gf = fs.createField<double>(
        atlas::option::name(varName)
        | atlas::option::levels(nLevels)
        | atlas::option::global());
    fs.gather(field, gf);
    globalFields.push_back(gf);
  }

  // --- 2. Root creates the NetCDF file, defines dims/vars, writes data ---
  if (isRoot) {
    std::remove(filepath.c_str());
    // NetCDF-4: classic formats cap a variable at 2-4 GB, and a 3D double
    // field reaches ~4.9 GB at 1/12 deg (41 x 3297 x 4500).
    int ncid;
    ncCheck(nc_create(filepath.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid), "nc_create", filepath);

    // Define dimensions
    int dim_time, dim_layer, dim_lath, dim_lonh;
    ncCheck(nc_def_dim(ncid, "Time",  NC_UNLIMITED, &dim_time), "nc_def_dim", filepath);
    ncCheck(nc_def_dim(ncid, "Layer", nz,           &dim_layer), "nc_def_dim", filepath);
    ncCheck(nc_def_dim(ncid, "lath",  nj,           &dim_lath), "nc_def_dim", filepath);
    ncCheck(nc_def_dim(ncid, "lonh",  ni,           &dim_lonh), "nc_def_dim", filepath);

    // Define variables (skip duplicates, e.g. two JEDI names → same file var)
    std::unordered_map<std::string, int> varids;
    fi = 0;
    for (const auto & field : x) {
      const std::string & varName = fileVarNames[fi++];
      if (varName.empty() || varids.count(varName)) continue;
      int varid;
      const int fIdx = fi - 1;
      const int dimY = nativeNorth.count(fIdx) ? dim_latq : dim_lath;
      const int dimX = nativePairs.count(fIdx) ? dim_lonq : dim_lonh;
      if (field.shape(1) > 1) {
        int dimids[4] = {dim_time, dim_layer, dim_lath, dim_lonh};
        ncCheck(nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 4, dimids, &varid),
                "nc_def_var " + varName, filepath);
      } else {
        int dimids[3] = {dim_time, dim_lath, dim_lonh};
        ncCheck(nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 3, dimids, &varid),
                "nc_def_var " + varName, filepath);
      }
      varids[varName] = varid;
    }
    ncCheck(nc_enddef(ncid), "nc_enddef", filepath);

    // Retrieve the gathered global indices (only valid on root)
    auto gIdx = atlas::array::make_view<double, 2>(globalGidxField);
    const int nGathered = static_cast<int>(globalGidxField.shape(0));

    // --- Diagnostic: validate gathered global indices ---
    {
      int nBadIdx = 0;
      double minIdx = 1.0e30, maxIdx = -1.0e30;
      for (int n = 0; n < nGathered; ++n) {
        const double g = gIdx(n, 0);
        if (g < minIdx) minIdx = g;
        if (g > maxIdx) maxIdx = g;
        const int flatIdx = static_cast<int>(g) - 1;
        if (flatIdx < 0 || flatIdx >= nj * ni) ++nBadIdx;
      }
      oops::Log::info() << "writeMOM6Netcdf: gathered " << nGathered
                        << " nodes (full grid = " << nj * ni << ")"
                        << ", gIdx range [" << minIdx << ", " << maxIdx << "]"
                        << ", out-of-range count = " << nBadIdx << std::endl;
    }

    // Expand a gathered field onto the full (nLevels, nj, ni) grid.
    auto expand = [&](const atlas::Field & gf) {
      auto gv = atlas::array::make_view<double, 2>(gf);
      const int nL = static_cast<int>(gf.shape(1));
      std::vector<double> full(static_cast<size_t>(nL) * nj * ni, 0.0);
      for (int n = 0; n < static_cast<int>(gf.shape(0)); ++n) {
        const size_t flatIdx = static_cast<size_t>(gIdx(n, 0)) - 1;
        for (int k = 0; k < nL; ++k)
          full[static_cast<size_t>(k) * nj * ni + flatIdx] = gv(n, k);
      }
      return full;
    };
    auto put = [&](const std::string & varName, int nL, const std::vector<double> & buf) {
      if (nL > 1) {
        const size_t start[4] = {0, 0, 0, 0};
        const size_t count[4] = {1, static_cast<size_t>(nL),
                                    static_cast<size_t>(nj),
                                    static_cast<size_t>(ni)};
        nc_put_vara_double(ncid, varids.at(varName), start, count, buf.data());
      } else {
        const size_t start[3] = {0, 0, 0};
        const size_t count[3] = {1, static_cast<size_t>(nj),
                                    static_cast<size_t>(ni)};
        nc_put_vara_double(ncid, varids.at(varName), start, count, buf.data());
      }
    };

    // Write each gathered field
    fi = 0;
    for (const auto & field : x) {
      const std::string & varName = fileVarNames[fi];
      const auto & gf = globalFields[fi];
      const int fIdx = fi;
      ++fi;
      if (varName.empty()) continue;

      if (nativeNorth.count(fIdx)) continue;  // written with its east partner
      if (nativePairs.count(fIdx)) {
        // Geographic at h-points -> grid-relative, then h-points -> C-grid faces
        const int iN = nativePairs.at(fIdx);
        const int nL = static_cast<int>(gf.shape(1));
        ASSERT(static_cast<int>(globalFields[iN].shape(1)) == nL);
        std::vector<double> uT = expand(gf);
        std::vector<double> vT = expand(globalFields[iN]);
        mom6::geoToGrid(velGrid->cosRot, velGrid->sinRot, nL, &uT, &vT);
        std::vector<double> u, v;
        mom6::tracerToFaces(uT, vT, velGrid->mask3d, ni, nj, nL, velGrid->hasFold, &u, &v);
        put(varName, nL, u);
        put(fileVarNames[iN], nL, v);
        oops::Log::info() << "writeMOM6Netcdf: " << field.name() << ", "
                          << x.field(iN).name() << " -> C-grid " << varName << ", "
                          << fileVarNames[iN] << " (rotated, interpolated to faces) to "
                          << filepath << std::endl;
        continue;
      }

      auto gView = atlas::array::make_view<double, 2>(gf);
      const int nGlobal = static_cast<int>(gf.shape(0));
      const int nLevels = static_cast<int>(gf.shape(1));

      // The gathered field has nActiveGlobal entries (ocean+fringe only).
      // Use the gathered global_index to place each value on the full grid.
      //   global_index = jG * ni + iG + 1  →  flatIdx = global_index - 1
      if (nLevels > 1) {
        const size_t nHoriz = static_cast<size_t>(nj) * ni;
        std::vector<double> buf(static_cast<size_t>(nz) * nHoriz, 0.0);
        for (int n = 0; n < nGlobal; ++n) {
          const size_t flatIdx = static_cast<size_t>(gIdx(n, 0)) - 1;
          for (int k = 0; k < nLevels; ++k)
            buf[static_cast<size_t>(k) * nHoriz + flatIdx] = gView(n, k);
        }
        const size_t start[4] = {0, 0, 0, 0};
        const size_t count[4] = {1, static_cast<size_t>(nz),
                                    static_cast<size_t>(nj),
                                    static_cast<size_t>(ni)};
        ncCheck(nc_put_vara_double(ncid, varids.at(varName), start, count, buf.data()),
                "nc_put_vara_double " + varName, filepath);
      } else {
        std::vector<double> buf(static_cast<size_t>(nj * ni), 0.0);
        for (int n = 0; n < nGlobal; ++n) {
          const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
          buf[flatIdx] = gView(n, 0);
        }
        const size_t start[3] = {0, 0, 0};
        const size_t count[3] = {1, static_cast<size_t>(nj),
                                    static_cast<size_t>(ni)};
        ncCheck(nc_put_vara_double(ncid, varids.at(varName), start, count, buf.data()),
                "nc_put_vara_double " + varName, filepath);
      }

      oops::Log::info() << "writeMOM6Netcdf: " << field.name()
                        << " -> " << varName << " to " << filepath << std::endl;
    }

    ncCheck(nc_close(ncid), "nc_close", filepath);
  }

  // Ensure the file is fully written before any rank returns (and potentially
  // tries to read it back).
  comm.barrier();
}

// ---------------------------------------------------------------------------
void writeMOM6NetcdfBlocks(const std::string & filepath,
                           const atlas::FieldSet & x,
                           const std::vector<std::string> & fileVarNames,
                           const int nz,
                           const GeometryMOM6 & geom) {
  const eckit::mpi::Comm & comm = geom.comm();
  const bool isRoot = (comm.rank() == 0);
  const mom6::Block block = geom.mom6Block();
  const size_t nLocal = block.size();

  if (isRoot) std::remove(filepath.c_str());
  comm.barrier();

  int ncid;
  ncCheck(nc_create_par(filepath.c_str(), NC_CLOBBER | NC_NETCDF4,
                        MPI_Comm_f2c(comm.communicator()), MPI_INFO_NULL, &ncid),
          "nc_create_par", filepath);

  int dim_time, dim_layer, dim_lath, dim_lonh;
  ncCheck(nc_def_dim(ncid, "Time",  NC_UNLIMITED, &dim_time), "nc_def_dim", filepath);
  ncCheck(nc_def_dim(ncid, "Layer", nz,           &dim_layer), "nc_def_dim", filepath);
  ncCheck(nc_def_dim(ncid, "lath",  geom.njEff(), &dim_lath), "nc_def_dim", filepath);
  ncCheck(nc_def_dim(ncid, "lonh",  geom.niEff(), &dim_lonh), "nc_def_dim", filepath);

  // Define variables (skip duplicates, e.g. two JEDI names → same file var)
  std::unordered_map<std::string, int> varids;
  int fi = 0;
  for (const auto & field : x) {
    const std::string & varName = fileVarNames[fi++];
    if (varName.empty() || varids.count(varName)) continue;
    int varid;
    if (field.shape(1) > 1) {
      int dimids[4] = {dim_time, dim_layer, dim_lath, dim_lonh};
      ncCheck(nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 4, dimids, &varid),
              "nc_def_var " + varName, filepath);
    } else {
      int dimids[3] = {dim_time, dim_lath, dim_lonh};
      ncCheck(nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 3, dimids, &varid),
              "nc_def_var " + varName, filepath);
    }
    varids[varName] = varid;
  }
  ncCheck(nc_enddef(ncid), "nc_enddef", filepath);

  // Each field: owned JEDI values → this task's MOM6 block (land stays 0) →
  // collective hyperslab write.
  fi = 0;
  for (const auto & field : x) {
    const std::string & varName = fileVarNames[fi++];
    if (varName.empty()) continue;
    const int nLevels = static_cast<int>(field.shape(1));

    atlas::Field mom6Field = geom.mom6FunctionSpace().createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(nLevels));
    auto vm = atlas::array::make_view<double, 2>(mom6Field);
    vm.assign(0.0);
    geom.gatherFromJedi(field, &mom6Field);

    std::vector<double> buf(static_cast<size_t>(nLevels) * nLocal);
    for (size_t n = 0; n < nLocal; ++n)
      for (int k = 0; k < nLevels; ++k)
        buf[static_cast<size_t>(k) * nLocal + n] = vm(n, k);

    const int varid = varids.at(varName);
    ncCheck(nc_var_par_access(ncid, varid, NC_COLLECTIVE), "nc_var_par_access", filepath);
    if (nLevels > 1) {
      const size_t start[4] = {0, 0, static_cast<size_t>(block.j0),
                               static_cast<size_t>(block.i0)};
      const size_t count[4] = {1, static_cast<size_t>(nLevels),
                               static_cast<size_t>(block.nj), static_cast<size_t>(block.ni)};
      ncCheck(nc_put_vara_double(ncid, varid, start, count, buf.data()),
              "nc_put_vara_double " + varName, filepath);
    } else {
      const size_t start[3] = {0, static_cast<size_t>(block.j0), static_cast<size_t>(block.i0)};
      const size_t count[3] = {1, static_cast<size_t>(block.nj), static_cast<size_t>(block.ni)};
      ncCheck(nc_put_vara_double(ncid, varid, start, count, buf.data()),
              "nc_put_vara_double " + varName, filepath);
    }
    if (isRoot)
      oops::Log::info() << "writeMOM6NetcdfBlocks: " << field.name()
                        << " -> " << varName << " to " << filepath << std::endl;
  }

  ncCheck(nc_close(ncid), "nc_close", filepath);
  comm.barrier();
}

// ---------------------------------------------------------------------------
}  // namespace ijedi

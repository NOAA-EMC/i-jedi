#include "ijedi/Io/mom6/ReadMOM6Netcdf.h"

#include <netcdf.h>

#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Logger.h"

namespace ijedi {

namespace {

std::string findDimName(int ncid,
                        const std::string & role,
                        const std::vector<std::string> & candidates,
                        size_t * len) {
  int dimid;
  for (const auto & name : candidates) {
    if (nc_inq_dimid(ncid, name.c_str(), &dimid) == NC_NOERR) {
      nc_inq_dimlen(ncid, dimid, len);
      return name;
    }
  }

  std::string tried;
  for (const auto & name : candidates) {
    if (!tried.empty()) tried += ", ";
    tried += name;
  }
  throw eckit::Exception(
      "readMOM6Netcdf: cannot find " + role +
          " dimension. Tried: " + tried,
      Here());
}

}  // namespace

// ---------------------------------------------------------------------------
void readMOM6Netcdf(const std::string & filepath,
                    atlas::FieldSet & x,
                    const std::vector<std::string> & fileVarNames,
                    const std::vector<double> & scalings,
                    const int numLevelsGeom,
                    const eckit::mpi::Comm & comm) {
  const bool isRoot = (comm.rank() == 0);
  const auto fs =
      atlas::functionspace::NodeColumns(x.field(0).functionspace());

  // --- 1. Gather global indices to root so we can map structured grid
  //        positions to the sparse global field. ---
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
  fs.gather(gidxField, globalGidxField);  // collective

  auto gIdx = atlas::array::make_view<double, 2>(globalGidxField);
  const int nGlobal = static_cast<int>(globalGidxField.shape(0));

  // --- 2. Root opens NetCDF and reads file dimensions. ---
  int ncid = -1;
  int file_ni = 0, file_nj = 0, file_nz = 0;
  if (isRoot) {
    ASSERT_MSG(nc_open(filepath.c_str(), NC_NOWRITE, &ncid) == NC_NOERR,
               "readMOM6Netcdf: cannot open " + filepath);

    size_t dimNi = 0, dimNj = 0, dimNz = 0;
    const std::string dimX =
      findDimName(ncid, "x (longitude)", {"lonh", "xh", "xaxis_1", "nx"}, &dimNi);
    const std::string dimY =
      findDimName(ncid, "y (latitude)", {"lath", "yh", "yaxis_1", "ny"}, &dimNj);
    const std::string dimZ =
      findDimName(ncid, "z (layer)", {"Layer", "z_l", "zaxis_1"}, &dimNz);

    file_ni = static_cast<int>(dimNi);
    file_nj = static_cast<int>(dimNj);
    file_nz = static_cast<int>(dimNz);

    if (file_nz != numLevelsGeom) {
      nc_close(ncid);
      throw eckit::BadValue(
        "readMOM6Netcdf: file dimensions in " + filepath +
          " are " + dimX + "=" + std::to_string(file_ni) +
          ", " + dimY + "=" + std::to_string(file_nj) +
          ", " + dimZ + "=" + std::to_string(file_nz) +
          " but geometry expects z=" + std::to_string(numLevelsGeom),
        Here());
    }

    const int nStructured = file_ni * file_nj;
    for (int n = 0; n < nGlobal; ++n) {
      const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
      if (flatIdx < 0 || flatIdx >= nStructured) {
        nc_close(ncid);
        throw eckit::BadValue(
            "readMOM6Netcdf: global index " + std::to_string(flatIdx + 1) +
                " is out of bounds for file horizontal dimensions " +
                dimX + "=" + std::to_string(file_ni) +
                ", " + dimY + "=" + std::to_string(file_nj) +
                " in " + filepath,
            Here());
      }
    }
  }

  // --- 3. For each field: root reads structured data into a global field,
  //        then scatter distributes to all ranks (collective). ---

  int fi = 0;
  for (auto & field : x) {
    const std::string & fileVarName = fileVarNames[fi];
    const double scale = scalings[fi];
    ++fi;

    const int nLevels = field.shape(1);

    // All ranks must create the global field (scatter is collective)
    atlas::Field gf = fs.createField<double>(
        atlas::option::name(field.name())
        | atlas::option::levels(nLevels)
        | atlas::option::global());

    if (isRoot && !fileVarName.empty()) {
      int varid;
      const bool found =
          (nc_inq_varid(ncid, fileVarName.c_str(), &varid) == NC_NOERR);

      if (found) {
        auto gView = atlas::array::make_view<double, 2>(gf);
        const bool is3D = (nLevels > 1);

        if (is3D) {
          std::vector<double> buf(
              static_cast<size_t>(file_nz) * file_nj * file_ni);
          const size_t start[4] = {0, 0, 0, 0};
          const size_t count[4] = {1, static_cast<size_t>(file_nz),
                                      static_cast<size_t>(file_nj),
                                      static_cast<size_t>(file_ni)};
          nc_get_vara_double(ncid, varid, start, count, buf.data());

          if (scale != 0.0)
            for (auto & v : buf) v *= scale;

          for (int n = 0; n < nGlobal; ++n) {
            const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
            const int iG = flatIdx % file_ni;
            const int jG = flatIdx / file_ni;
            for (int k = 0; k < nLevels; ++k) {
              gView(n, k) = buf[static_cast<size_t>(k) * file_nj * file_ni
                              + static_cast<size_t>(jG) * file_ni
                              + static_cast<size_t>(iG)];
            }
          }
        } else {
          std::vector<double> buf(
              static_cast<size_t>(file_nj) * file_ni);
          const size_t start[3] = {0, 0, 0};
          const size_t count[3] = {1, static_cast<size_t>(file_nj),
                                      static_cast<size_t>(file_ni)};
          nc_get_vara_double(ncid, varid, start, count, buf.data());

          if (scale != 0.0)
            for (auto & v : buf) v *= scale;

          for (int n = 0; n < nGlobal; ++n) {
            const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
            gView(n, 0) = buf[flatIdx];
          }
        }

        oops::Log::info() << "readMOM6Netcdf: " << field.name()
                          << " <- " << fileVarName << " from " << filepath
                          << std::endl;
      } else {
        oops::Log::warning() << "readMOM6Netcdf: variable '" << fileVarName
                             << "' not found in file, skipping '" << field.name()
                             << "'" << std::endl;
      }
    }

    fs.scatter(gf, field);  // collective
  }

  if (isRoot) nc_close(ncid);

  // Halo exchange to fill ghost nodes from their owning ranks
  for (auto & field : x) {
    atlas::functionspace::NodeColumns(field.functionspace()).haloExchange(field);
  }
}

// ---------------------------------------------------------------------------
}  // namespace ijedi

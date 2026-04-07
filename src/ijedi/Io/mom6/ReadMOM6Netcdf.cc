#include "ijedi/Io/mom6/ReadMOM6Netcdf.h"

#include <netcdf.h>

#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace ijedi {

// ---------------------------------------------------------------------------
void readMOM6Netcdf(const std::string & filepath,
                    atlas::FieldSet & x,
                    const std::vector<std::string> & fileVarNames,
                    const std::vector<double> & scalings,
                    const int numLevelsGeom) {
  // Open netCDF file (read-only, all ranks read independently)
  int ncid;
  ASSERT_MSG(nc_open(filepath.c_str(), NC_NOWRITE, &ncid) == NC_NOERR,
             "readMOM6Netcdf: cannot open " + filepath);

  // Read horizontal and vertical file dimensions.
  // Convention: global_index = jG * lonh + iG + 1  (1-based)
  int dimid;
  size_t file_ni, file_nj, file_nz;
  ASSERT_MSG(nc_inq_dimid(ncid, "lonh",  &dimid) == NC_NOERR,
             "readMOM6Netcdf: missing dimension 'lonh' in " + filepath);
  nc_inq_dimlen(ncid, dimid, &file_ni);
  ASSERT_MSG(nc_inq_dimid(ncid, "lath",  &dimid) == NC_NOERR,
             "readMOM6Netcdf: missing dimension 'lath' in " + filepath);
  nc_inq_dimlen(ncid, dimid, &file_nj);
  ASSERT_MSG(nc_inq_dimid(ncid, "Layer", &dimid) == NC_NOERR,
             "readMOM6Netcdf: missing dimension 'Layer' in " + filepath);
  nc_inq_dimlen(ncid, dimid, &file_nz);

  // Sanity-check vertical dimension against geometry
  ASSERT_MSG(static_cast<int>(file_nz) == numLevelsGeom,
             "readMOM6Netcdf: file has " + std::to_string(file_nz)
             + " layers but geometry has " + std::to_string(numLevelsGeom));

  // Global index view for scattering structured data into unstructured nodes
  const auto fs =
      atlas::functionspace::NodeColumns(x.field(0).functionspace());
  auto g_view = atlas::array::make_view<atlas::gidx_t, 1>(
      fs.nodes().global_index());
  const int nNodes = static_cast<int>(fs.nodes().size());

  int fi = 0;
  for (auto & field : x) {
    const std::string & fileVarName = fileVarNames[fi];
    const double scale = scalings[fi];
    ++fi;

    if (fileVarName.empty()) continue;

    int varid;
    if (nc_inq_varid(ncid, fileVarName.c_str(), &varid) != NC_NOERR) {
      oops::Log::warning() << "readMOM6Netcdf: variable '" << fileVarName
                           << "' not found in file, skipping '" << field.name()
                           << "'" << std::endl;
      continue;
    }

    auto view = atlas::array::make_view<double, 2>(field);
    const int nLevels = field.shape(1);
    const bool is3D = (nLevels > 1);

    if (is3D) {
      // Read full (1, nz, nj, ni) slab at Time=0
      std::vector<double> buf(file_nz * file_nj * file_ni);
      const size_t start[4] = {0, 0, 0, 0};
      const size_t count[4] = {1, file_nz, file_nj, file_ni};
      nc_get_vara_double(ncid, varid, start, count, buf.data());

      if (scale != 0.0)
        for (auto & v : buf) v *= scale;

      for (int n = 0; n < nNodes; ++n) {
        const int flatIdx = static_cast<int>(g_view(n) - 1);
        const int iG = flatIdx % static_cast<int>(file_ni);
        const int jG = flatIdx / static_cast<int>(file_ni);
        for (int k = 0; k < nLevels; ++k) {
          view(n, k) = buf[static_cast<size_t>(k) * file_nj * file_ni
                         + static_cast<size_t>(jG) * file_ni
                         + static_cast<size_t>(iG)];
        }
      }
    } else {
      // 2D variable: read full (1, nj, ni) slab at Time=0
      std::vector<double> buf(file_nj * file_ni);
      const size_t start[3] = {0, 0, 0};
      const size_t count[3] = {1, file_nj, file_ni};
      nc_get_vara_double(ncid, varid, start, count, buf.data());

      if (scale != 0.0)
        for (auto & v : buf) v *= scale;

      for (int n = 0; n < nNodes; ++n) {
        const int flatIdx = static_cast<int>(g_view(n) - 1);
        view(n, 0) = buf[flatIdx];
      }
    }

    oops::Log::info() << "readMOM6Netcdf: " << field.name()
                      << " <- " << fileVarName << " from " << filepath << std::endl;
  }

  nc_close(ncid);

  // Halo exchange to fill ghost nodes from their owning ranks
  for (auto & field : x) {
    atlas::functionspace::NodeColumns(field.functionspace()).haloExchange(field);
  }
}

// ---------------------------------------------------------------------------
}  // namespace ijedi

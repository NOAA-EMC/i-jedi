#include "ijedi/Io/mpas/ReadMPASNetcdf.h"

#include <netcdf.h>

#include <numeric>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace ijedi {
namespace {

std::string dimName(const int ncid, const int dimid) {
  char name[NC_MAX_NAME + 1];
  size_t len = 0;
  const int status = nc_inq_dim(ncid, dimid, name, &len);
  ASSERT_MSG(status == NC_NOERR, std::string("readMPASNetcdf: nc_inq_dim failed: ") +
                                 nc_strerror(status));
  return std::string(name);
}

std::vector<size_t> dimLengths(const int ncid, const std::vector<int> & dimids) {
  std::vector<size_t> lengths(dimids.size(), 0);
  for (size_t i = 0; i < dimids.size(); ++i) {
    const int status = nc_inq_dimlen(ncid, dimids[i], &lengths[i]);
    ASSERT_MSG(status == NC_NOERR, std::string("readMPASNetcdf: nc_inq_dimlen failed: ") +
                                   nc_strerror(status));
  }
  return lengths;
}

std::vector<size_t> makeStrides(const std::vector<size_t> & dims) {
  std::vector<size_t> strides(dims.size(), 1);
  if (dims.empty()) {
    return strides;
  }
  for (int i = static_cast<int>(dims.size()) - 2; i >= 0; --i) {
    strides[i] = strides[i + 1] * dims[i + 1];
  }
  return strides;
}

size_t linearIndex(const std::vector<size_t> & indices, const std::vector<size_t> & strides) {
  size_t idx = 0;
  for (size_t i = 0; i < indices.size(); ++i) {
    idx += indices[i] * strides[i];
  }
  return idx;
}

}  // namespace

void readMPASNetcdf(const std::string & filepath,
                    atlas::FieldSet & x,
                    const std::vector<std::string> & fileVarNames,
                    const std::vector<double> & scalings,
                    const int nCellsGlobal,
                    const eckit::mpi::Comm & comm) {
  const bool isRoot = (comm.rank() == 0);
  const auto fs = atlas::functionspace::NodeColumns(x.field(0).functionspace());

  auto localGidx = atlas::array::make_view<atlas::gidx_t, 1>(fs.nodes().global_index());
  atlas::Field gidxField = fs.createField<double>(atlas::option::name("__gidx") |
                                                   atlas::option::levels(1));
  {
    auto gv = atlas::array::make_view<double, 2>(gidxField);
    for (atlas::idx_t n = 0; n < fs.nodes().size(); ++n) {
      gv(n, 0) = static_cast<double>(localGidx(n));
    }
  }

  atlas::Field globalGidxField = fs.createField<double>(atlas::option::name("__gidx") |
                                                         atlas::option::levels(1) |
                                                         atlas::option::global());
  fs.gather(gidxField, globalGidxField);

  int ncid = -1;
  if (isRoot) {
    const int status = nc_open(filepath.c_str(), NC_NOWRITE, &ncid);
    ASSERT_MSG(status == NC_NOERR,
               std::string("readMPASNetcdf: cannot open ") + filepath + ": " +
                   nc_strerror(status));
  }

  auto gIdx = atlas::array::make_view<double, 2>(globalGidxField);
  const int nGlobal = static_cast<int>(globalGidxField.shape(0));

  int fi = 0;
  for (auto & field : x) {
    const std::string & fileVarName = fileVarNames.at(fi);
    const double scale = scalings.at(fi);
    ++fi;

    const int nLevelsField = field.shape(1);

    atlas::Field gf = fs.createField<double>(atlas::option::name(field.name()) |
                                             atlas::option::levels(nLevelsField) |
                                             atlas::option::global());

    if (isRoot && !fileVarName.empty()) {
      int varid = -1;
      const int found = nc_inq_varid(ncid, fileVarName.c_str(), &varid);
      if (found == NC_NOERR) {
        char varName[NC_MAX_NAME + 1];
        nc_type xtype = NC_NAT;
        int ndims = 0;
        int dimids[NC_MAX_VAR_DIMS] = {0};
        int natts = 0;

        const int inqStatus = nc_inq_var(ncid, varid, varName, &xtype, &ndims, dimids, &natts);
        ASSERT_MSG(inqStatus == NC_NOERR,
                   std::string("readMPASNetcdf: nc_inq_var failed for ") + fileVarName + ": " +
                       nc_strerror(inqStatus));

        std::vector<int> dimidVec(dimids, dimids + ndims);
        std::vector<size_t> dims = dimLengths(ncid, dimidVec);

        int cellPos = -1;
        int levelPos = -1;
        int timePos = -1;

        for (int i = 0; i < ndims; ++i) {
          const std::string dn = dimName(ncid, dimids[i]);
          if (dn == "nCells") {
            cellPos = i;
          } else if (dn == "nVertLevels" || dn == "nVertLevelsP1") {
            levelPos = i;
          } else if (dn == "Time") {
            timePos = i;
          }
        }

        ASSERT_MSG(cellPos >= 0,
                   "readMPASNetcdf: variable '" + fileVarName + "' has no nCells dimension");
        ASSERT_MSG(static_cast<int>(dims[cellPos]) == nCellsGlobal,
                   "readMPASNetcdf: variable '" + fileVarName +
                       "' nCells does not match geometry nCellsGlobal");

        if (nLevelsField > 1) {
          ASSERT_MSG(levelPos >= 0,
                     "readMPASNetcdf: variable '" + fileVarName +
                         "' has no vertical dimension for a 3D JEDI field");
          ASSERT_MSG(static_cast<int>(dims[levelPos]) == nLevelsField,
                     "readMPASNetcdf: variable '" + fileVarName +
                         "' vertical size does not match field levels");
        }

        for (int i = 0; i < ndims; ++i) {
          if (i == cellPos || i == levelPos || i == timePos) {
            continue;
          }
          ASSERT_MSG(dims[i] == 1,
                     "readMPASNetcdf: unsupported extra dimension with length > 1 in variable '" +
                         fileVarName + "'");
        }

        size_t total = 1;
        for (const size_t d : dims) {
          total *= d;
        }

        std::vector<double> buffer(total, 0.0);
        const int getStatus = nc_get_var_double(ncid, varid, buffer.data());
        ASSERT_MSG(getStatus == NC_NOERR,
                   std::string("readMPASNetcdf: nc_get_var_double failed for ") + fileVarName +
                       ": " + nc_strerror(getStatus));

        if (scale != 0.0) {
          for (double & v : buffer) {
            v *= scale;
          }
        }

        const std::vector<size_t> strides = makeStrides(dims);
        auto gView = atlas::array::make_view<double, 2>(gf);

        for (int n = 0; n < nGlobal; ++n) {
          const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
          ASSERT_MSG(flatIdx >= 0 && flatIdx < nCellsGlobal,
                     "readMPASNetcdf: gathered global index out of range");

          if (nLevelsField > 1) {
            for (int k = 0; k < nLevelsField; ++k) {
              std::vector<size_t> index(dims.size(), 0);
              index[cellPos] = static_cast<size_t>(flatIdx);
              if (timePos >= 0) {
                index[timePos] = 0;
              }
              index[levelPos] = static_cast<size_t>(k);
              gView(n, k) = buffer[linearIndex(index, strides)];
            }
          } else {
            std::vector<size_t> index(dims.size(), 0);
            index[cellPos] = static_cast<size_t>(flatIdx);
            if (timePos >= 0) {
              index[timePos] = 0;
            }
            if (levelPos >= 0) {
              index[levelPos] = 0;
            }
            gView(n, 0) = buffer[linearIndex(index, strides)];
          }
        }

        oops::Log::info() << "readMPASNetcdf: " << field.name() << " <- " << fileVarName
                          << " from " << filepath << std::endl;
      } else {
        oops::Log::warning() << "readMPASNetcdf: variable '" << fileVarName
                             << "' not found in file, skipping '" << field.name() << "'"
                             << std::endl;
      }
    }

    fs.scatter(gf, field);
    field.metadata().set("interp_type", "default");
  }

  if (isRoot) {
    nc_close(ncid);
  }

  for (auto & field : x) {
    atlas::functionspace::NodeColumns(field.functionspace()).haloExchange(field);
  }
}

}  // namespace ijedi

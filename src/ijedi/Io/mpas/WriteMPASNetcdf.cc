#include "ijedi/Io/mpas/WriteMPASNetcdf.h"

#include <netcdf.h>

#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
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
  ASSERT_MSG(status == NC_NOERR, std::string("writeMPASNetcdf: nc_inq_dim failed: ") +
                                 nc_strerror(status));
  return std::string(name);
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

int ensureDimension(const int ncid, const std::string & name, const size_t len) {
  int dimid = -1;
  const int inq = nc_inq_dimid(ncid, name.c_str(), &dimid);
  if (inq == NC_NOERR) {
    size_t existingLen = 0;
    const int lenStatus = nc_inq_dimlen(ncid, dimid, &existingLen);
    ASSERT_MSG(lenStatus == NC_NOERR,
               std::string("writeMPASNetcdf: nc_inq_dimlen failed for ") + name + ": " +
                   nc_strerror(lenStatus));
    ASSERT_MSG(existingLen == len,
               "writeMPASNetcdf: existing dim '" + name + "' has incompatible length");
    return dimid;
  }

  const int defStatus = nc_def_dim(ncid, name.c_str(), len, &dimid);
  ASSERT_MSG(defStatus == NC_NOERR,
             std::string("writeMPASNetcdf: nc_def_dim failed for ") + name + ": " +
                 nc_strerror(defStatus));
  return dimid;
}

int ensureVariable(const int ncid,
                   const std::string & varName,
                   const int nLevels,
                   const int nVertLevels,
                   const int dimTime,
                   const int dimCells,
                   const int dimVert,
                   const int dimVertP1) {
  int varid = -1;
  if (nc_inq_varid(ncid, varName.c_str(), &varid) == NC_NOERR) {
    return varid;
  }

  if (nLevels > 1) {
    int levelDim = dimVert;
    if (nLevels == nVertLevels + 1) {
      levelDim = dimVertP1;
    }
    int dimids[3] = {dimTime, dimCells, levelDim};
    const int status = nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 3, dimids, &varid);
    ASSERT_MSG(status == NC_NOERR,
               std::string("writeMPASNetcdf: nc_def_var failed for ") + varName + ": " +
                   nc_strerror(status));
  } else {
    int dimids[2] = {dimTime, dimCells};
    const int status = nc_def_var(ncid, varName.c_str(), NC_DOUBLE, 2, dimids, &varid);
    ASSERT_MSG(status == NC_NOERR,
               std::string("writeMPASNetcdf: nc_def_var failed for ") + varName + ": " +
                   nc_strerror(status));
  }

  return varid;
}

}  // namespace

void writeMPASNetcdf(const std::string & filepath,
                     const atlas::FieldSet & x,
                     const std::vector<std::string> & fileVarNames,
                     const int nCellsGlobal,
                     const int nVertLevels,
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

  std::vector<atlas::Field> globalFields;
  globalFields.reserve(x.size());

  int fi = 0;
  for (const auto & field : x) {
    const std::string & varName = fileVarNames.at(fi++);
    if (varName.empty()) {
      globalFields.emplace_back();
      continue;
    }
    const int nLevels = field.shape(1);
    atlas::Field gf = fs.createField<double>(atlas::option::name(varName) |
                                             atlas::option::levels(nLevels) |
                                             atlas::option::global());
    fs.gather(field, gf);
    globalFields.push_back(gf);
  }

  if (isRoot) {
    int ncid = -1;
    int openStatus = nc_open(filepath.c_str(), NC_WRITE, &ncid);
    if (openStatus != NC_NOERR) {
      std::remove(filepath.c_str());
      openStatus = nc_create(filepath.c_str(), NC_CLOBBER, &ncid);
      ASSERT_MSG(openStatus == NC_NOERR,
                 std::string("writeMPASNetcdf: cannot create ") + filepath + ": " +
                     nc_strerror(openStatus));
    }

    int dimTime = -1;
    int dimCells = -1;
    int dimVert = -1;
    int dimVertP1 = -1;

    const int redefStatus = nc_redef(ncid);
    ASSERT_MSG(redefStatus == NC_NOERR || redefStatus == NC_EINDEFINE,
               std::string("writeMPASNetcdf: nc_redef failed: ") + nc_strerror(redefStatus));

    if (nc_inq_dimid(ncid, "Time", &dimTime) != NC_NOERR) {
      const int status = nc_def_dim(ncid, "Time", NC_UNLIMITED, &dimTime);
      ASSERT_MSG(status == NC_NOERR,
                 std::string("writeMPASNetcdf: nc_def_dim failed for Time: ") +
                     nc_strerror(status));
    }

    dimCells = ensureDimension(ncid, "nCells", static_cast<size_t>(nCellsGlobal));
    dimVert = ensureDimension(ncid, "nVertLevels", static_cast<size_t>(nVertLevels));
    dimVertP1 = ensureDimension(ncid, "nVertLevelsP1", static_cast<size_t>(nVertLevels + 1));

    std::unordered_map<std::string, int> varids;
    fi = 0;
    for (const auto & field : x) {
      const std::string & varName = fileVarNames.at(fi++);
      if (varName.empty() || varids.count(varName) > 0) {
        continue;
      }
      varids[varName] = ensureVariable(ncid, varName, field.shape(1), nVertLevels,
                                       dimTime, dimCells, dimVert, dimVertP1);
    }

    const int enddefStatus = nc_enddef(ncid);
    ASSERT_MSG(enddefStatus == NC_NOERR,
               std::string("writeMPASNetcdf: nc_enddef failed: ") + nc_strerror(enddefStatus));

    auto gIdx = atlas::array::make_view<double, 2>(globalGidxField);

    fi = 0;
    for (const auto & field : x) {
      const std::string & varName = fileVarNames.at(fi);
      const auto & gf = globalFields.at(fi);
      ++fi;

      if (varName.empty()) {
        continue;
      }

      int varid = varids.at(varName);

      char varNameInFile[NC_MAX_NAME + 1];
      nc_type xtype = NC_NAT;
      int ndims = 0;
      int dimids[NC_MAX_VAR_DIMS] = {0};
      int natts = 0;
      const int inqVarStatus = nc_inq_var(ncid, varid, varNameInFile, &xtype, &ndims, dimids, &natts);
      ASSERT_MSG(inqVarStatus == NC_NOERR,
                 std::string("writeMPASNetcdf: nc_inq_var failed for ") + varName + ": " +
                     nc_strerror(inqVarStatus));

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
                 "writeMPASNetcdf: variable '" + varName + "' has no nCells dimension");

      const int nLevelsField = field.shape(1);
      if (nLevelsField > 1) {
        ASSERT_MSG(levelPos >= 0,
                   "writeMPASNetcdf: variable '" + varName +
                       "' has no vertical dimension for a 3D JEDI field");
        size_t levelDimLen = 0;
        nc_inq_dimlen(ncid, dimids[levelPos], &levelDimLen);
        ASSERT_MSG(static_cast<int>(levelDimLen) == nLevelsField,
                   "writeMPASNetcdf: variable '" + varName +
                       "' vertical size does not match field levels");
      }

      // Build start/count for nc_put_vara — always write exactly one time record (index 0).
      // This avoids relying on nc_inq_dimlen for the UNLIMITED Time dim, which returns 0
      // on a freshly created file (causing a zero-size buffer and heap overflow otherwise).
      std::vector<size_t> start(ndims, 0);
      std::vector<size_t> count(ndims, 1);
      for (int i = 0; i < ndims; ++i) {
        if (i == cellPos) {
          count[i] = static_cast<size_t>(nCellsGlobal);
        } else if (i == levelPos) {
          count[i] = static_cast<size_t>(nLevelsField);
        }
        // timePos stays start=0, count=1
      }

      // Allocate a compact (nCells × nLevels) buffer in C-order: [cell][level]
      const size_t bufTotal = static_cast<size_t>(nCellsGlobal) *
                              static_cast<size_t>(nLevelsField > 1 ? nLevelsField : 1);
      std::vector<double> buffer(bufTotal, 0.0);

      auto gView = atlas::array::make_view<double, 2>(gf);
      const int nGlobal = static_cast<int>(gf.shape(0));

      // Build strides for the compact buffer (same dimension ordering as count, minus Time).
      // The compact layout is always [nCells][nLevels] regardless of the on-file dim order.
      // We'll rearrange into on-file order when filling the buffer.
      // Actually: match the on-file layout so nc_put_vara writes directly.
      // On-file order for start/count is the same as the dim order of the variable.
      // We need to fill the buffer in that same order.
      std::vector<size_t> bufDims(ndims);
      for (int i = 0; i < ndims; ++i) {
        bufDims[i] = count[i];
      }
      const std::vector<size_t> bufStrides = makeStrides(bufDims);

      for (int n = 0; n < nGlobal; ++n) {
        const int flatIdx = static_cast<int>(gIdx(n, 0)) - 1;
        ASSERT_MSG(flatIdx >= 0 && flatIdx < nCellsGlobal,
                   "writeMPASNetcdf: gathered global index out of range");

        if (nLevelsField > 1) {
          for (int k = 0; k < nLevelsField; ++k) {
            std::vector<size_t> index(ndims, 0);
            index[cellPos] = static_cast<size_t>(flatIdx);
            index[levelPos] = static_cast<size_t>(k);
            // timePos index stays 0
            buffer[linearIndex(index, bufStrides)] = gView(n, k);
          }
        } else {
          std::vector<size_t> index(ndims, 0);
          index[cellPos] = static_cast<size_t>(flatIdx);
          buffer[linearIndex(index, bufStrides)] = gView(n, 0);
        }
      }

      const int putStatus = nc_put_vara_double(ncid, varid, start.data(), count.data(),
                                               buffer.data());
      ASSERT_MSG(putStatus == NC_NOERR,
                 std::string("writeMPASNetcdf: nc_put_vara_double failed for ") + varName +
                     ": " + nc_strerror(putStatus));

      oops::Log::info() << "writeMPASNetcdf: " << field.name() << " -> " << varName
                        << " to " << filepath << std::endl;
    }

    nc_close(ncid);
  }

  comm.barrier();
}

}  // namespace ijedi

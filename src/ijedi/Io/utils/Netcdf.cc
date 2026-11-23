#include "ijedi/Io/utils/Netcdf.h"

#include <netcdf.h>

#include <algorithm>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace ijedi {

namespace {

// Where the horizontal position and the levels of a field are in a variable: for each layout
// dimension, its length and its stride in the C-order variable buffer
struct VariableMap {
  int varid = -1;
  size_t size = 0;                                    // number of values of the variable
  std::vector<size_t> hLen, hStride, lLen, lStride;
};

// Offset in the variable buffer of the flat index `flat` over dimensions (len, stride), outer to
// inner
size_t offset(size_t flat, const std::vector<size_t> & len, const std::vector<size_t> & stride) {
  size_t off = 0;
  for (size_t d = len.size(); d-- > 0;) {
    off += (flat % len[d]) * stride[d];
    flat /= len[d];
  }
  return off;
}

// Rank 0: maps a variable onto the layout, or returns an error message
std::string mapVariable(const int ncid, const std::string & varName, const int nLevels,
                        const NetcdfLayout & layout, const size_t maxFlat, VariableMap & map) {
  const std::string what = "variable " + varName + ": ";
  if (nc_inq_varid(ncid, varName.c_str(), &map.varid) != NC_NOERR) return "no " + what;
  int ndims = 0;
  std::vector<int> dimids(NC_MAX_VAR_DIMS);
  if (nc_inq_varndims(ncid, map.varid, &ndims) != NC_NOERR
      || nc_inq_vardimid(ncid, map.varid, dimids.data()) != NC_NOERR) {
    return what + "cannot inquire its dimensions";
  }
  std::vector<std::string> names(ndims);
  std::vector<size_t> lens(ndims), strides(ndims);
  for (int d = 0; d < ndims; ++d) {
    char name[NC_MAX_NAME + 1];
    if (nc_inq_dim(ncid, dimids[d], name, &lens[d]) != NC_NOERR)
      return what + "cannot inquire its dimensions";
    names[d] = name;
  }
  map.size = 1;
  for (int d = ndims; d-- > 0;) {
    strides[d] = map.size;
    map.size *= lens[d];
  }
  auto find = [&](const std::string & dim) {
    return std::find(names.begin(), names.end(), dim) - names.begin();
  };
  size_t nHoriz = 1, nLev = 1;
  for (const auto & dim : layout.horizontalDims) {
    const auto d = find(dim);
    if (d == ndims) return what + "no " + dim + " dimension";
    map.hLen.push_back(lens[d]);
    map.hStride.push_back(strides[d]);
    nHoriz *= lens[d];
  }
  for (const auto & dim : layout.levelDims) {
    const auto d = find(dim);
    if (d == ndims) {
      if (nLevels == 1) continue;   // a single-level field may come from a variable without it
      return what + "no " + dim + " dimension";
    }
    map.lLen.push_back(lens[d]);
    map.lStride.push_back(strides[d]);
    nLev *= lens[d];
  }
  for (int d = 0; d < ndims; ++d) {
    const bool inLayout =
      std::count(layout.horizontalDims.begin(), layout.horizontalDims.end(), names[d])
      + std::count(layout.levelDims.begin(), layout.levelDims.end(), names[d]) > 0;
    if (!inLayout && lens[d] != 1)
      return what + "dimension " + names[d] + " has length " + std::to_string(lens[d])
             + " (dimensions outside the layout must have length 1)";
  }
  if (nLev != static_cast<size_t>(nLevels))
    return what + std::to_string(nLev) + " levels in the file, " + std::to_string(nLevels)
           + " in the field";
  if (maxFlat >= nHoriz)
    return what + "global index " + std::to_string(maxFlat + 1) + " beyond the "
           + std::to_string(nHoriz) + " horizontal points of the file";
  return "";
}

// The gathered global indices (rank 0) and the variables of an open file, mapped onto the layout
struct OpenFile {
  atlas::Field gidx;            // global field: global index of every gathered node
  int nGathered = 0;            // number of gathered nodes (rank 0), 0 elsewhere
  int ncid = -1;                // rank 0
  std::vector<VariableMap> maps;
};

// Gathers the global indices, then rank 0 opens the file and maps every variable before any
// further collective call. A failure is broadcast so that all ranks throw together: a throw on
// rank 0 alone would leave the other ranks waiting in the collectives that follow.
OpenFile openFile(const std::string & caller, const std::string & filepath, const int mode,
                  const atlas::FieldSet & x, const std::vector<std::string> & fileVarNames,
                  const NetcdfLayout & layout, const eckit::mpi::Comm & comm) {
  ASSERT(fileVarNames.size() == static_cast<size_t>(x.size()));
  const bool isRoot = (comm.rank() == 0);
  const atlas::FunctionSpace fs = x.field(0).functionspace();
  OpenFile file;
  file.maps.resize(fileVarNames.size());

  const auto localGidx = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());
  atlas::Field gidxField = fs.createField<double>(
      atlas::option::name("__gidx") | atlas::option::levels(1));
  {
    auto gv = atlas::array::make_view<double, 2>(gidxField);
    for (atlas::idx_t n = 0; n < gidxField.shape(0); ++n) gv(n, 0) = localGidx(n);
  }
  file.gidx = fs.createField<double>(
      atlas::option::name("__gidx") | atlas::option::levels(1) | atlas::option::global());
  fs.gather(gidxField, file.gidx);
  const auto gIdx = atlas::array::make_view<double, 2>(file.gidx);
  file.nGathered = isRoot ? static_cast<int>(file.gidx.shape(0)) : 0;

  std::string error;
  if (isRoot) {
    size_t maxFlat = 0;
    for (int n = 0; n < file.nGathered; ++n) {
      if (gIdx(n, 0) < 1) {
        error = "invalid global index " + std::to_string(gIdx(n, 0));
        break;
      }
      maxFlat = std::max(maxFlat, static_cast<size_t>(gIdx(n, 0)) - 1);
    }
    if (error.empty() && nc_open(filepath.c_str(), mode, &file.ncid) != NC_NOERR) {
      file.ncid = -1;
      error = "cannot open " + filepath
              + (mode == NC_WRITE ? " (it is updated in place and must exist)" : "");
    }
    for (size_t f = 0; f < fileVarNames.size() && error.empty(); ++f) {
      if (fileVarNames[f].empty()) continue;
      error = mapVariable(file.ncid, fileVarNames[f], x.field(f).shape(1), layout, maxFlat,
                          file.maps[f]);
      if (!error.empty()) error += " in " + filepath;
    }
    if (!error.empty() && file.ncid >= 0) nc_close(file.ncid);
  }
  int failed = error.empty() ? 0 : 1;
  comm.broadcast(failed, 0);
  if (failed) {
    throw eckit::UserError(caller + ": "
                           + (isRoot ? error : std::string("check failed on rank 0")), Here());
  }
  return file;
}

// An i/o error past the checks of openFile aborts all ranks (rank 0 cannot throw alone in the
// middle of the collectives)
void abortOnError(const int status, const std::string & what, const eckit::mpi::Comm & comm) {
  if (status != NC_NOERR) {
    oops::Log::error() << what << ": " << nc_strerror(status) << std::endl;
    comm.abort(1);
  }
}

// Offsets of the levels of a field in the variable buffer
std::vector<size_t> levelOffsets(const VariableMap & map, const int nLevels) {
  std::vector<size_t> levelOffset(nLevels, 0);
  if (!map.lLen.empty()) {
    for (int k = 0; k < nLevels; ++k) levelOffset[k] = offset(k, map.lLen, map.lStride);
  }
  return levelOffset;
}

}  // namespace

// -------------------------------------------------------------------------------------------------

void updateNetcdf(const std::string & filepath,
                  const atlas::FieldSet & x,
                  const std::vector<std::string> & fileVarNames,
                  const NetcdfLayout & layout,
                  const eckit::mpi::Comm & comm) {
  if (x.size() == 0) return;
  const bool isRoot = (comm.rank() == 0);
  const atlas::FunctionSpace fs = x.field(0).functionspace();
  const OpenFile file = openFile("updateNetcdf", filepath, NC_WRITE, x, fileVarNames, layout,
                                 comm);
  const auto gIdx = atlas::array::make_view<double, 2>(file.gidx);

  // One field at a time: gather it on root (collective), which updates the variable at the
  // gathered nodes, so that root holds a single global field at a time. A failed read never
  // writes a stale buffer over the values JEDI does not hold.
  for (size_t f = 0; f < fileVarNames.size(); ++f) {
    if (fileVarNames[f].empty()) continue;
    const atlas::Field & field = x.field(f);
    const int nLevels = field.shape(1);
    atlas::Field gf = fs.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(nLevels)
        | atlas::option::global());
    fs.gather(field, gf);
    if (!isRoot) continue;

    const VariableMap & map = file.maps[f];
    std::vector<double> buf(map.size);
    abortOnError(nc_get_var_double(file.ncid, map.varid, buf.data()),
                 "updateNetcdf: cannot read " + fileVarNames[f] + " in " + filepath, comm);
    const std::vector<size_t> levelOffset = levelOffsets(map, nLevels);
    const auto gView = atlas::array::make_view<double, 2>(gf);
    for (int n = 0; n < file.nGathered; ++n) {
      const size_t h = offset(static_cast<size_t>(gIdx(n, 0)) - 1, map.hLen, map.hStride);
      for (int k = 0; k < nLevels; ++k) buf[h + levelOffset[k]] = gView(n, k);
    }
    abortOnError(nc_put_var_double(file.ncid, map.varid, buf.data()),
                 "updateNetcdf: cannot write " + fileVarNames[f] + " in " + filepath, comm);
    oops::Log::info() << "updateNetcdf: " << field.name() << " -> " << fileVarNames[f]
                      << " in " << filepath << std::endl;
  }
  if (isRoot) abortOnError(nc_close(file.ncid), "updateNetcdf: cannot close " + filepath, comm);
  comm.barrier();
}

// -------------------------------------------------------------------------------------------------

void readNetcdf(const std::string & filepath,
                atlas::FieldSet & x,
                const std::vector<std::string> & fileVarNames,
                const NetcdfLayout & layout,
                const eckit::mpi::Comm & comm) {
  if (x.size() == 0) return;
  const bool isRoot = (comm.rank() == 0);
  const atlas::FunctionSpace fs = x.field(0).functionspace();
  const OpenFile file = openFile("readNetcdf", filepath, NC_NOWRITE, x, fileVarNames, layout,
                                 comm);
  const auto gIdx = atlas::array::make_view<double, 2>(file.gidx);

  // One field at a time: root fills a global field in the gather order, which the scatter
  // (collective) distributes; the halo exchange then sets the ghost nodes
  for (size_t f = 0; f < fileVarNames.size(); ++f) {
    if (fileVarNames[f].empty()) continue;
    atlas::Field & field = x.field(f);
    const int nLevels = field.shape(1);
    atlas::Field gf = fs.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(nLevels)
        | atlas::option::global());
    if (isRoot) {
      const VariableMap & map = file.maps[f];
      std::vector<double> buf(map.size);
      abortOnError(nc_get_var_double(file.ncid, map.varid, buf.data()),
                   "readNetcdf: cannot read " + fileVarNames[f] + " in " + filepath, comm);
      double fill;
      if (nc_get_att_double(file.ncid, map.varid, "_FillValue", &fill) == NC_NOERR) {
        for (auto & v : buf) if (v == fill) v = util::missingValue<double>();
      }
      const std::vector<size_t> levelOffset = levelOffsets(map, nLevels);
      auto gView = atlas::array::make_view<double, 2>(gf);
      for (int n = 0; n < file.nGathered; ++n) {
        const size_t h = offset(static_cast<size_t>(gIdx(n, 0)) - 1, map.hLen, map.hStride);
        for (int k = 0; k < nLevels; ++k) gView(n, k) = buf[h + levelOffset[k]];
      }
    }
    fs.scatter(gf, field);
    fs.haloExchange(field);
    oops::Log::info() << "readNetcdf: " << fileVarNames[f] << " in " << filepath << " -> "
                      << field.name() << std::endl;
  }
  if (isRoot) abortOnError(nc_close(file.ncid), "readNetcdf: cannot close " + filepath, comm);
}

}  // namespace ijedi

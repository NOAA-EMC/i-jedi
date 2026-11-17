#include "ijedi/Io/cice/UpdateCICENetcdf.h"

#include <netcdf.h>

#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace ijedi {

void updateCICENetcdf(const std::string & filepath,
                      const atlas::FieldSet & x,
                      const std::vector<std::string> & fileVarNames,
                      const int ni, const int nj,
                      const eckit::mpi::Comm & comm) {
  const bool isRoot = (comm.rank() == 0);
  const auto fs = atlas::functionspace::NodeColumns(x.field(0).functionspace());

  // Global index (j * ni + i + 1) of every gathered node
  auto localGidx = atlas::array::make_view<atlas::gidx_t, 1>(fs.nodes().global_index());
  atlas::Field gidxField = fs.createField<double>(
      atlas::option::name("__gidx") | atlas::option::levels(1));
  {
    auto gv = atlas::array::make_view<double, 2>(gidxField);
    for (atlas::idx_t n = 0; n < fs.nodes().size(); ++n) gv(n, 0) = localGidx(n);
  }
  atlas::Field globalGidxField = fs.createField<double>(
      atlas::option::name("__gidx") | atlas::option::levels(1) | atlas::option::global());
  fs.gather(gidxField, globalGidxField);

  const size_t nHoriz = static_cast<size_t>(ni) * nj;

  // Rank 0 checks the file and every variable before any collective call. A failure is
  // broadcast so that all ranks throw together: a throw on rank 0 alone would leave the other
  // ranks waiting in the gathers below.
  int ncid = -1;
  std::vector<int> varids(fileVarNames.size(), -1);
  std::string error;
  if (isRoot) {
    if (nc_open(filepath.c_str(), NC_WRITE, &ncid) != NC_NOERR) {
      error = "cannot open " + filepath
              + " (the restart is updated in place and must be staged beforehand)";
    }
    int fi = 0;
    for (const auto & field : x) {
      const std::string & varName = fileVarNames[fi++];
      if (!error.empty()) break;
      if (varName.empty()) continue;
      int varid, ndims;
      if (nc_inq_varid(ncid, varName.c_str(), &varid) != NC_NOERR) {
        error = "no variable " + varName + " in " + filepath;
        break;
      }
      std::vector<int> dimids(NC_MAX_VAR_DIMS);
      if (nc_inq_varndims(ncid, varid, &ndims) != NC_NOERR
          || nc_inq_vardimid(ncid, varid, dimids.data()) != NC_NOERR) {
        error = "cannot inquire " + varName + " in " + filepath;
        break;
      }
      size_t nTotal = 1;
      std::vector<size_t> lens(ndims);
      for (int d = 0; d < ndims && error.empty(); ++d) {
        if (nc_inq_dimlen(ncid, dimids[d], &lens[d]) != NC_NOERR) {
          error = "cannot inquire " + varName + " in " + filepath;
        }
        nTotal *= lens[d];
      }
      if (!error.empty()) break;
      const int nLevels = field.shape(1);
      if (ndims < 2 || lens[ndims - 1] != static_cast<size_t>(ni)
          || lens[ndims - 2] != static_cast<size_t>(nj) || nTotal != nHoriz * nLevels) {
        error = "shape of " + varName + " in " + filepath + " does not match a field with "
                + std::to_string(nLevels) + " levels";
        break;
      }
      varids[fi - 1] = varid;
    }
    if (!error.empty() && ncid >= 0) nc_close(ncid);
  }
  int failed = error.empty() ? 0 : 1;
  comm.broadcast(failed, 0);
  if (failed) {
    throw eckit::UserError("updateCICENetcdf: "
                           + (isRoot ? error : std::string("check failed on rank 0")),
                           Here());
  }

  const auto gIdx = atlas::array::make_view<double, 2>(globalGidxField);
  const int nGathered = static_cast<int>(globalGidxField.shape(0));

  // One field at a time: gather it on root (collective), which overwrites the file variable at
  // the mesh nodes, so that root holds a single global field at a time. An i/o error past the
  // checks above aborts all ranks (rank 0 cannot throw alone in the middle of the gathers), and
  // a failed read never writes a stale buffer over the values JEDI does not hold.
  auto abortOnError = [&](const int status, const std::string & what) {
    if (status != NC_NOERR) {
      oops::Log::error() << "updateCICENetcdf: " << what << ": " << nc_strerror(status)
                         << std::endl;
      comm.abort(1);
    }
  };
  int fi = 0;
  for (const auto & field : x) {
    const std::string & varName = fileVarNames[fi++];
    if (varName.empty()) continue;
    atlas::Field gf = fs.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(field.shape(1))
        | atlas::option::global());
    fs.gather(field, gf);
    if (!isRoot) continue;

    const int varid = varids[fi - 1];
    const int nLevels = field.shape(1);
    std::vector<double> buf(nHoriz * nLevels);
    abortOnError(nc_get_var_double(ncid, varid, buf.data()),
                 "cannot read " + varName + " in " + filepath);
    const auto gView = atlas::array::make_view<double, 2>(gf);
    for (int n = 0; n < nGathered; ++n) {
      const size_t flat = static_cast<size_t>(gIdx(n, 0)) - 1;
      ASSERT(flat < nHoriz);
      for (int k = 0; k < nLevels; ++k) buf[k * nHoriz + flat] = gView(n, k);
    }
    abortOnError(nc_put_var_double(ncid, varid, buf.data()),
                 "cannot write " + varName + " in " + filepath);
    oops::Log::info() << "updateCICENetcdf: " << field.name() << " -> " << varName
                      << " in " << filepath << std::endl;
  }
  if (isRoot) abortOnError(nc_close(ncid), "cannot close " + filepath);
  comm.barrier();
}

}  // namespace ijedi

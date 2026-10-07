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

  // Gather the fields (collective)
  std::vector<atlas::Field> globalFields;
  int fi = 0;
  for (const auto & field : x) {
    if (fileVarNames[fi++].empty()) {
      globalFields.emplace_back();
      continue;
    }
    atlas::Field gf = fs.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(field.shape(1))
        | atlas::option::global());
    fs.gather(field, gf);
    globalFields.push_back(gf);
  }

  if (isRoot) {
    int ncid;
    ASSERT_MSG(nc_open(filepath.c_str(), NC_WRITE, &ncid) == NC_NOERR,
               "updateCICENetcdf: cannot open " + filepath);
    const auto gIdx = atlas::array::make_view<double, 2>(globalGidxField);
    const int nGathered = static_cast<int>(globalGidxField.shape(0));
    const size_t nHoriz = static_cast<size_t>(ni) * nj;

    fi = 0;
    for (const auto & field : x) {
      const std::string & varName = fileVarNames[fi];
      const atlas::Field & gf = globalFields[fi++];
      if (varName.empty()) continue;

      int varid;
      ASSERT_MSG(nc_inq_varid(ncid, varName.c_str(), &varid) == NC_NOERR,
                 "updateCICENetcdf: no variable " + varName + " in " + filepath);
      int ndims;
      nc_inq_varndims(ncid, varid, &ndims);
      std::vector<int> dimids(ndims);
      nc_inq_vardimid(ncid, varid, dimids.data());
      size_t nTotal = 1;
      std::vector<size_t> lens(ndims);
      for (int d = 0; d < ndims; ++d) {
        nc_inq_dimlen(ncid, dimids[d], &lens[d]);
        nTotal *= lens[d];
      }
      const int nLevels = field.shape(1);
      ASSERT_MSG(ndims >= 2 && lens[ndims - 1] == static_cast<size_t>(ni)
                 && lens[ndims - 2] == static_cast<size_t>(nj)
                 && nTotal == nHoriz * nLevels,
                 "updateCICENetcdf: shape of " + varName + " in " + filepath
                 + " does not match a field with " + std::to_string(nLevels) + " levels");

      // Read the variable, overwrite the mesh nodes, write it back
      std::vector<double> buf(nTotal);
      nc_get_var_double(ncid, varid, buf.data());
      const auto gView = atlas::array::make_view<double, 2>(gf);
      for (int n = 0; n < nGathered; ++n) {
        const size_t flat = static_cast<size_t>(gIdx(n, 0)) - 1;
        ASSERT(flat < nHoriz);
        for (int k = 0; k < nLevels; ++k) buf[k * nHoriz + flat] = gView(n, k);
      }
      ASSERT_MSG(nc_put_var_double(ncid, varid, buf.data()) == NC_NOERR,
                 "updateCICENetcdf: cannot write " + varName + " in " + filepath);
      oops::Log::info() << "updateCICENetcdf: " << field.name() << " -> " << varName
                        << " in " << filepath << std::endl;
    }
    nc_close(ncid);
  }
  comm.barrier();
}

}  // namespace ijedi

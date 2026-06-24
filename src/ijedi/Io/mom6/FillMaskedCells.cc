#include "ijedi/Io/mom6/FillMaskedCells.h"

#include <queue>
#include <string>
#include <unordered_set>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/HybridElements.h"

#include "oops/util/Logger.h"

namespace ijedi {

void applyBoundaryConditions(atlas::FieldSet & x,
                             const atlas::Field & mask2d,
                             const FieldsMetadata & meta) {
  const auto & longNames = meta.getLongNames();
  std::unordered_set<std::string> tracerNames;
  for (const auto & name : longNames) {
    if (meta.getFieldMetadata(name).getIsTracer()) tracerNames.insert(name);
  }
  const std::unordered_set<std::string> knownNames(longNames.begin(), longNames.end());

  const auto fs  = atlas::functionspace::NodeColumns(mask2d.functionspace());
  const int npts = static_cast<int>(mask2d.shape(0));
  auto maskView  = atlas::array::make_view<double, 2>(mask2d);

  // Build node-to-node adjacency from cell-node connectivity (once, reused per field).
  const auto & conn   = fs.mesh().cells().node_connectivity();
  const int    nCells = static_cast<int>(conn.rows());
  std::vector<std::vector<int>> adj(npts);
  for (int c = 0; c < nCells; ++c) {
    const int nv = static_cast<int>(conn.cols(c));
    for (int a = 0; a < nv; ++a) {
      for (int b = a + 1; b < nv; ++b) {
        const int na = static_cast<int>(conn(c, a));
        const int nb = static_cast<int>(conn(c, b));
        if (na < npts && nb < npts) {
          adj[na].push_back(nb);
          adj[nb].push_back(na);
        }
      }
    }
  }

  for (auto & field : x) {
    if (field.shape(1) != 1) continue;  // 3-D fields deferred

    const std::string & fname = field.name();
    if (!knownNames.count(fname)) {
      oops::Log::warning() << "applyBoundaryConditions: '" << fname
                           << "' not in metadata, skipping" << std::endl;
      continue;
    }

    auto fView = atlas::array::make_view<double, 2>(field);

    if (!tracerNames.count(fname)) {
      // Non-tracer: no-flux / no-slip → zero all masked nodes.
      int nZeroed = 0;
      for (int n = 0; n < npts; ++n) {
        if (maskView(n, 0) <= 0.5) { fView(n, 0) = 0.0; ++nZeroed; }
      }
      oops::Log::info() << "applyBoundaryConditions: " << fname
                        << " zero BC on " << nZeroed << " masked nodes" << std::endl;
      continue;
    }

    // Tracer: extrapolate each masked (land) node from its nearest ocean node,
    // imposing a Neumann / zero-normal-gradient BC. Implemented as a flood fill
    // outward from every ocean node at once, so each land node takes the value
    // of the closest ocean node by mesh connectivity.
    //
    // Ghost ocean nodes are included as seeds — they already carry valid values
    // from the scatter + halo exchange in readMOM6Netcdf, so coast nodes adjacent
    // to another rank's ocean domain are seeded without extra communication.
    std::vector<bool> filled(npts, false);
    std::queue<int>   q;
    for (int n = 0; n < npts; ++n) {
      if (maskView(n, 0) > 0.5) { filled[n] = true; q.push(n); }
    }

    int nFilled = 0;
    while (!q.empty()) {
      const int cur = q.front(); q.pop();
      for (const int nb : adj[cur]) {
        if (!filled[nb]) {
          filled[nb]   = true;
          fView(nb, 0) = fView(cur, 0);
          q.push(nb);
          ++nFilled;
        }
      }
    }

    // Propagate flooded owned-node values into ghost copies on neighbouring ranks.
    fs.haloExchange(field);

    oops::Log::info() << "applyBoundaryConditions: " << fname
                      << " Neumann (flood fill) on " << nFilled
                      << " masked nodes" << std::endl;
  }
}

}  // namespace ijedi

#include "ijedi/Io/mom6/FillMaskedCells.h"

#include <queue>
#include <string>
#include <unordered_set>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/option.h"

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
  // All node pairs sharing a cell are considered adjacent.
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

  // Auxiliary field to carry fill-status across MPI ranks via halo exchange.
  atlas::Field filledFld = fs.createField<double>(
      atlas::option::name("__bfs_filled") | atlas::option::levels(1));
  auto filledView = atlas::array::make_view<double, 2>(filledFld);

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

    // Tracer: BFS flood fill from ocean nodes (Neumann / zero normal gradient BC).
    // Nodes seeded in index order → deterministic on each rank.
    for (int n = 0; n < npts; ++n)
      filledView(n, 0) = (maskView(n, 0) > 0.5) ? 1.0 : 0.0;

    std::vector<bool> localFilled(npts, false);
    std::queue<int>   q;
    for (int n = 0; n < npts; ++n) {
      if (maskView(n, 0) > 0.5) { localFilled[n] = true; q.push(n); }
    }

    // Alternate local BFS and halo exchange until convergence.
    // Each outer iteration pushes the filled frontier one halo-width further.
    int nHalo = 0;
    for (;;) {
      while (!q.empty()) {
        const int cur = q.front(); q.pop();
        for (const int nb : adj[cur]) {
          if (!localFilled[nb]) {
            localFilled[nb]   = true;
            filledView(nb, 0) = 1.0;
            fView(nb, 0)      = fView(cur, 0);
            q.push(nb);
          }
        }
      }

      // Propagate filled values and fill-status to ghost nodes on other ranks.
      fs.haloExchange(field);
      fs.haloExchange(filledFld);
      ++nHalo;

      // Reseed from ghost nodes whose owners have now filled them.
      bool newSeeds = false;
      for (int n = 0; n < npts; ++n) {
        if (!localFilled[n] && filledView(n, 0) > 0.5) {
          localFilled[n] = true;
          q.push(n);
          newSeeds = true;
        }
      }
      if (!newSeeds) break;
    }

    oops::Log::info() << "applyBoundaryConditions: " << fname
                      << " Neumann (BFS flood) after " << nHalo
                      << " halo exchange(s)" << std::endl;
  }
}

}  // namespace ijedi

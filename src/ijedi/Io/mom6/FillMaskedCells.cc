#include "ijedi/Io/mom6/FillMaskedCells.h"

#include <algorithm>
#include <array>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "atlas/array.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace ijedi {

void applyBoundaryConditions(std::vector<atlas::Field> & fields,
                             const atlas::Field & mask3d,
                             const std::vector<atlas::gidx_t> & gidx,
                             const int ni, const int nj, const bool hasFold,
                             const FieldsMetadata & meta) {
  const auto & longNames = meta.getLongNames();
  const std::unordered_set<std::string> knownNames(longNames.begin(), longNames.end());

  const int npts  = static_cast<int>(mask3d.shape(0));
  const int nmask = static_cast<int>(mask3d.shape(1));
  ASSERT(static_cast<int>(gidx.size()) == npts);
  auto maskView   = atlas::array::make_view<const double, 2>(mask3d);

  // Rows in global index order, and the row of each global node
  std::vector<int> order(npts);
  for (int n = 0; n < npts; ++n) order[n] = n;
  std::sort(order.begin(), order.end(), [&](int a, int b) { return gidx[a] < gidx[b]; });
  std::unordered_map<atlas::gidx_t, int> rowOf;
  for (int n = 0; n < npts; ++n) rowOf[gidx[n]] = n;
  const auto row = [&](int i, int j) {
    const auto it = rowOf.find(static_cast<atlas::gidx_t>(j) * ni + i + 1);
    return it == rowOf.end() ? -1 : it->second;
  };

  // Cells of the geometry's mesh (GeometryMOM6): quads whose four corners are mesh nodes,
  // periodic in i, then the tripolar fold quads on the last row
  std::vector<std::array<int, 4>> cells;
  for (const int n : order) {
    const int i = static_cast<int>((gidx[n] - 1) % ni);
    const int j = static_cast<int>((gidx[n] - 1) / ni);
    if (j + 1 >= nj) continue;
    const int i1 = (i + 1) % ni;
    const std::array<int, 4> c = {n, row(i1, j), row(i1, j + 1), row(i, j + 1)};
    if (std::min({c[1], c[2], c[3]}) >= 0) cells.push_back(c);
  }
  if (hasFold) {
    for (const int n : order) {
      const int i = static_cast<int>((gidx[n] - 1) % ni);
      const int j = static_cast<int>((gidx[n] - 1) / ni);
      if (j != nj - 1 || 2 * (i + 1) >= ni || ni - 2 - i < 0) continue;
      const std::array<int, 4> c = {n, row(i + 1, j), row(ni - 2 - i, j), row(ni - 1 - i, j)};
      if (std::min({c[1], c[2], c[3]}) >= 0) cells.push_back(c);
    }
  }

  // Node-to-node adjacency: every pair of corners of a cell
  std::vector<std::vector<int>> adj(npts);
  for (const auto & c : cells) {
    for (int a = 0; a < 4; ++a) {
      for (int b = a + 1; b < 4; ++b) {
        adj[c[a]].push_back(c[b]);
        adj[c[b]].push_back(c[a]);
      }
    }
  }

  for (auto & field : fields) {
    const std::string & fname = field.name();
    if (!knownNames.count(fname)) {
      oops::Log::warning() << "applyBoundaryConditions: '" << fname
                           << "' not in metadata, skipping" << std::endl;
      continue;
    }

    const std::string bcType = meta.getFieldMetadata(fname).getBcType();

    // "none": leave the field untouched (e.g. geometry-provided coordinates).
    if (bcType == "none") {
      oops::Log::info() << "applyBoundaryConditions: " << fname
                        << " bctype=none, no boundary condition applied" << std::endl;
      continue;
    }

    const int nlev  = static_cast<int>(field.shape(1));
    auto fView = atlas::array::make_view<double, 2>(field);

    // For each level k, the mask value at (n, k) drives the BC — this covers
    // both 2-D fields (nlev == 1, k always 0) and 3-D fields uniformly.
    // mask3d has nmask levels; clamp k to the last mask level so that a field
    // with more levels than the mask (shouldn't happen in practice) is safe.

    if (bcType == "zero") {
      // "zero": no-flux / no-slip → zero every masked node at every level.
      int nZeroed = 0;
      for (int n = 0; n < npts; ++n) {
        for (int k = 0; k < nlev; ++k) {
          const int km = std::min(k, nmask - 1);
          if (maskView(n, km) <= 0.5) { fView(n, k) = 0.0; ++nZeroed; }
        }
      }
      oops::Log::info() << "applyBoundaryConditions: " << fname
                        << " zero BC on " << nZeroed << " masked (node,level) pairs"
                        << " (" << nlev << " levels)" << std::endl;
      continue;
    }

    // "extrapolate": flood-fill each level independently from ocean seeds at that
    // level (Neumann / zero-gradient boundary condition).
    // Levels are processed top-down so that any node the horizontal flood cannot
    // reach (a mesh component with no wet node at this level, e.g. below-bottom
    // layers) can fall back to the value directly above it, which is guaranteed
    // filled by induction. At the surface (k == 0) there is no level above, so
    // an unreachable node is set to 0 as a last resort.
    int nFilled = 0;
    int nVertical = 0;
    for (int k = 0; k < nlev; ++k) {
      const int km = std::min(k, nmask - 1);
      std::vector<bool> filled(npts, false);
      std::queue<int>   q;
      for (const int n : order) {
        if (maskView(n, km) > 0.5) { filled[n] = true; q.push(n); }
      }
      while (!q.empty()) {
        const int cur = q.front(); q.pop();
        for (const int nb : adj[cur]) {
          if (!filled[nb]) {
            filled[nb]   = true;
            fView(nb, k) = fView(cur, k);
            q.push(nb);
            ++nFilled;
          }
        }
      }
      // Vertical fallback for nodes the horizontal flood never reached.
      for (int n = 0; n < npts; ++n) {
        if (!filled[n]) {
          fView(n, k) = (k > 0) ? fView(n, k - 1) : 0.0;
          ++nVertical;
        }
      }
    }

    oops::Log::info() << "applyBoundaryConditions: " << fname
                      << " Neumann (flood fill) on " << nFilled
                      << " masked nodes, vertical fallback on " << nVertical
                      << " masked nodes (" << nlev << " levels)" << std::endl;
  }
}

}  // namespace ijedi

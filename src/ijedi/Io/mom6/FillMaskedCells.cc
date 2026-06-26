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
                             const atlas::Field & mask3d,
                             const FieldsMetadata & meta) {
  const auto & longNames = meta.getLongNames();
  const std::unordered_set<std::string> knownNames(longNames.begin(), longNames.end());

  const auto fs   = atlas::functionspace::NodeColumns(mask3d.functionspace());
  const int npts  = static_cast<int>(mask3d.shape(0));
  const int nmask = static_cast<int>(mask3d.shape(1));
  auto maskView   = atlas::array::make_view<double, 2>(mask3d);

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
    // Ghost ocean nodes are included as seeds — they carry valid values from the
    // scatter + halo exchange in readMOM6Netcdf.
    //
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
      for (int n = 0; n < npts; ++n) {
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

    // Propagate flooded owned-node values into ghost copies on neighbouring ranks.
    fs.haloExchange(field);

    oops::Log::info() << "applyBoundaryConditions: " << fname
                      << " Neumann (flood fill) on " << nFilled
                      << " masked nodes (" << nlev << " levels)" << std::endl;
  }
}

}  // namespace ijedi

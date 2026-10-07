#pragma once

#include <vector>

#include "atlas/field/Field.h"
#include "atlas/library/config.h"

#include "ijedi/FieldMetadata/FieldsMetadata.h"

namespace ijedi {

/// Apply boundary conditions to masked (land) cells after reading.
///
/// The boundary condition applied to each field is controlled by the per-field
/// `bctype` attribute in the FieldsMetadata registry:
///
///   - "extrapolate": masked cells are filled by nearest-neighbour flood-fill
///     from the closest unmasked ocean node at the same level, approximating a
///     zero normal gradient (Neumann) boundary condition.  Levels are processed
///     top-down; any node unreachable by horizontal flood falls back to the
///     value at the level above (or 0 at the surface).
///
///   - "zero": masked cells are set to zero, imposing a no-flux / no-slip
///     boundary condition.
///
///   - "none": the field is left untouched (e.g. geometry-provided coordinates).
///
/// Both 2-D (nLevels == 1) and 3-D fields are processed uniformly; fields
/// absent from the metadata registry are skipped.
///
/// The fill works on fields gathered on one task (one row per global node, as
/// readMOM6Netcdf holds them before the scatter). The node graph is the one of
/// the geometry's mesh cells, rebuilt from the global indices: quads whose four
/// corners are mesh nodes (periodic in i) plus, on tripolar grids, the fold
/// quads along the last row. Seeds, cells and neighbours are taken in global
/// index order, so the filled values do not depend on the MPI decomposition.
///
/// @param fields   Global fields (rows = global nodes), filled in place.
/// @param mask3d   Global 3-D mask (1 = ocean, 0 = masked/land), same rows.
/// @param gidx     1-based global index (j * ni + i + 1) of each row.
/// @param ni, nj   Global grid dimensions.
/// @param hasFold  True for a tripolar grid with a northern fold.
/// @param meta     Field metadata registry providing the per-field bctype.
void applyBoundaryConditions(std::vector<atlas::Field> & fields,
                             const atlas::Field & mask3d,
                             const std::vector<atlas::gidx_t> & gidx,
                             const int ni, const int nj, const bool hasFold,
                             const FieldsMetadata & meta);

}  // namespace ijedi

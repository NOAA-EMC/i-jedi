#pragma once

#include "atlas/field/FieldSet.h"

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
/// The node adjacency graph is built once from the mesh connectivity and
/// reused across all fields. The per-level mask from mask3d drives the fill.
///
/// @param x      FieldSet to fill in-place (NodeColumns function space).
/// @param mask3d Geometry 3-D mask field (1 = ocean, 0 = masked/land).
/// @param meta   Field metadata registry providing the per-field bctype.
void applyBoundaryConditions(atlas::FieldSet & x,
                             const atlas::Field & mask3d,
                             const FieldsMetadata & meta);

}  // namespace ijedi

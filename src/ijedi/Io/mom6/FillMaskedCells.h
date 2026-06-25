#pragma once

#include "atlas/field/FieldSet.h"

#include "ijedi/FieldMetadata/FieldsMetadata.h"

namespace ijedi {

/// Apply boundary conditions to masked (land) cells after reading.
///
/// For tracer fields (isTracer == true in FieldsMetadata), masked cells are
/// filled by nearest-neighbour extrapolation from the closest unmasked ocean
/// node, approximating a zero normal gradient (Neumann) boundary condition.
///
/// For non-tracer fields (e.g. velocity components), masked cells are set to
/// zero, imposing a no-flux / no-slip boundary condition.
///
/// Only 2-D surface fields (nLevels == 1) are processed; 3-D fields and fields
/// absent from the metadata registry are skipped.
///
/// The KD-tree of ocean nodes is built once from mask2d and reused across all
/// tracer fields in the set.
///
/// @param x      FieldSet to fill in-place (NodeColumns function space).
/// @param mask2d Geometry mask field (1 = ocean, 0 = masked/land).
/// @param meta   Field metadata registry for tracer classification.
void applyBoundaryConditions(atlas::FieldSet & x,
                             const atlas::Field & mask2d,
                             const FieldsMetadata & meta);

}  // namespace ijedi

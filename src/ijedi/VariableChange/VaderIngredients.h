#pragma once

namespace atlas {
class FieldSet;
}  // namespace atlas

namespace ijedi {

class Geometry;

/// \brief Inject geometry-sourced ingredient fields required by several Vader
///        ocean recipes (e.g. SeaWaterTemperature_A/_B) that are not state
///        variables.
///
/// \details The fields are cloned from the (already built) geometry, so no grid
///          is regenerated. They are added only if not already present:
///            - latitude          <- geometry "lat"  (single level)
///            - longitude         <- geometry "lon"  (single level)
///            - sea_area_fraction <- geometry "mask2d", broadcast across all
///              model levels (SeaWaterTemperature_B indexes it per level)
void addVaderGeometryIngredients(atlas::FieldSet & fset, const Geometry & geom);

}  // namespace ijedi

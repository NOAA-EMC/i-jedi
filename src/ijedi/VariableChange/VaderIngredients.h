#pragma once

namespace atlas {
class FieldSet;
}  // namespace atlas

namespace ijedi {

class Geometry;

/// \brief Inject geometry-sourced ingredient fields required by several Vader
///        recipes (e.g. SeaWaterTemperature_A/_B) that are not state variables.
///
/// \details Added only if not already present, and split by responsibility:
///            - latitude / longitude: derived generically from the function
///              space (atlas lonlat()), so this works for any model regardless
///              of how its geometry names coordinate fields.
///            - model-specific ingredients (e.g. sea_area_fraction for the
///              ocean): delegated to Geometry::addModelVaderIngredients, keeping
///              model-specific logic out of this shared code.
void addVaderGeometryIngredients(atlas::FieldSet & fset, const Geometry & geom);

}  // namespace ijedi

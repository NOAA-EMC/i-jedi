#include "ijedi/VariableChange/VaderIngredients.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/option.h"

#include "ijedi/Geometry/Geometry.h"

namespace ijedi {

void addVaderGeometryIngredients(atlas::FieldSet & fset, const Geometry & geom) {
  // latitude / longitude: sourced generically from the function space so this
  // works for any model regardless of how its geometry names coordinate fields.
  // atlas lonlat() packs [lon, lat] into a single 2-component field.
  const atlas::FunctionSpace & fs = geom.functionSpace();
  if (!fset.has("latitude") || !fset.has("longitude")) {
    const auto lonlat = atlas::array::make_view<double, 2>(fs.lonlat());
    atlas::Field lonF = fs.createField<double>(
        atlas::option::name("longitude") | atlas::option::levels(1));
    atlas::Field latF = fs.createField<double>(
        atlas::option::name("latitude") | atlas::option::levels(1));
    auto lonView = atlas::array::make_view<double, 2>(lonF);
    auto latView = atlas::array::make_view<double, 2>(latF);
    for (atlas::idx_t n = 0; n < lonF.shape(0); ++n) {
      lonView(n, 0) = lonlat(n, 0);
      latView(n, 0) = lonlat(n, 1);
    }
    if (!fset.has("longitude")) fset.add(lonF);
    if (!fset.has("latitude"))  fset.add(latF);
  }

  // Model-specific ingredients (masks, area fractions, ...) are provided by the
  // geometry implementation so ocean-specific logic stays out of shared code.
  geom.addModelVaderIngredients(fset);
}

}  // namespace ijedi

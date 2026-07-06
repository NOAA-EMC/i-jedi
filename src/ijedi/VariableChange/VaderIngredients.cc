#include "ijedi/VariableChange/VaderIngredients.h"

#include <string>
#include <utility>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/option.h"

#include "ijedi/Geometry/Geometry.h"

namespace ijedi {

void addVaderGeometryIngredients(atlas::FieldSet & fset, const Geometry & geom) {
  // latitude / longitude: single-level coordinates cloned directly from the
  // cached geometry fields.
  const std::pair<std::string, std::string> coordMap[] = {
      {"latitude", "lat"}, {"longitude", "lon"}};
  for (const auto & [vaderName, geomName] : coordMap) {
    if (!fset.has(vaderName) && geom.fields().has(geomName)) {
      atlas::Field f = geom.fields().field(geomName).clone();
      f.rename(vaderName);
      fset.add(f);
    }
  }

  // sea_area_fraction: required by SeaWaterTemperature_B's Jacobian, which
  // indexes it per level. The geometry only carries a 2D land/sea mask, so
  // broadcast it across all model levels to keep the column access in bounds.
  if (!fset.has("sea_area_fraction") && geom.fields().has("mask2d")) {
    const atlas::Field & mask2d = geom.fields().field("mask2d");
    const int nlevels = geom.numLevels();
    atlas::Field saf = mask2d.functionspace().createField<double>(
        atlas::option::name("sea_area_fraction") | atlas::option::levels(nlevels));
    const auto maskView = atlas::array::make_view<double, 2>(mask2d);
    auto safView = atlas::array::make_view<double, 2>(saf);
    for (atlas::idx_t n = 0; n < saf.shape(0); ++n) {
      for (int k = 0; k < nlevels; ++k) {
        safView(n, k) = maskView(n, 0);
      }
    }
    fset.add(saf);
  }
}

}  // namespace ijedi

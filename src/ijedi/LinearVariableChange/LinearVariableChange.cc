/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/LinearVariableChange/LinearVariableChange.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/State.h"
#include "ijedi/VariableChange/VaderCookbook.h"
#include "ijedi/VariableChange/VaderIngredients.h"
#include "mist/base/ModelData.h"
#include "oops/base/Variables.h"

namespace ijedi {

LinearVariableChange::LinearVariableChange(const Geometry & geometry,
                                           const eckit::Configuration & varchangeConfig)
    : geom_(geometry) {
  eckit::LocalConfiguration configCookbook{};
  // If config has "vader cookbook" then use that, else use default
  if (varchangeConfig.has("vader cookbook")) {
    configCookbook = varchangeConfig.getSubConfiguration("vader cookbook");
  } else {
    const auto cb = ijedi::vaderDefaultCookbook();
    for (const auto & [key, val] : cb) {
      configCookbook.set(key, val);
    }
  }
  eckit::LocalConfiguration config{};
  config.set(vader::configCookbookKey, configCookbook);
  // Set up model data for cookbook
  auto configModelData = mist::base::ModelData(geometry).modelData();
  config.set(vader::configModelVarsKey, configModelData);

  initVaderVariableChange(configCookbook, configModelData);
}

void LinearVariableChange::changeVarTraj(const State & xx, const oops::Variables & vars) {
  // SeaWaterTemperature_B's Jacobian needs latitude, longitude and
  // sea_area_fraction in the trajectory fieldset (vader then computes the
  // trajectory sea_water_temperature itself via the NL recipe). These are
  // geometry coordinates/masks, not state variables, so inject them into a
  // working copy of the trajectory before setting the linearization point.
  State traj(xx);
  addVaderGeometryIngredients(traj.fieldSet(), geom_);
  mist::base::LinearVariableChange::changeVarTraj(traj, vars);
}

}  // namespace ijedi

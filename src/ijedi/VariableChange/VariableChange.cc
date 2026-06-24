/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/VariableChange/VariableChange.h"

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "ijedi/VariableChange/VaderCookbook.h"
#include "ijedi/VariableChange/VaderIngredients.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/State.h"
#include "mist/base/ModelData.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "vader/vader.h"

namespace ijedi {

VariableChange::VariableChange(const eckit::Configuration & varchangeConfig,
                               const Geometry & geometry) : geom_(geometry) {
  oops::Log::trace() << "ijedi::VariableChange::VariableChange starting" << std::endl;
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

  varchange_ = std::make_unique<mist::utils::VariableChange>(config);
  oops::Log::trace() << "ijedi::VariableChange::VariableChange done" << std::endl;
}

void VariableChange::changeVar(State & xx, const oops::Variables & vars) const {
  oops::Log::trace() << "ijedi::VariableChange::changeVar starting" << std::endl;

  // Several Vader recipes (e.g. SeaWaterTemperature_A) require geometry-sourced
  // ingredient fields (latitude, longitude, sea_area_fraction) that are not
  // state variables. Inject them transiently into the working FieldSet;
  // mist::utils::VariableChange filters its output back down to the requested
  // variables, so they do not persist in xx. Skip injection when no transform
  // is needed (mist returns early in that case, which would otherwise leave the
  // fields in the state).
  if (!(vars == xx.variables())) {
    addVaderGeometryIngredients(xx.fieldSet(), geom_);
  }

  varchange_->changeVar(xx, vars);
  xx.setAtlasFieldMetadata();
  oops::Log::trace() << "ijedi::VariableChange::changeVar done" << std::endl;
}

void VariableChange::changeVarInverse(State & xx, const oops::Variables & vars) const {
  throw eckit::NotImplemented("ijedi::VariableChange::changeVarInverse is not implemented", Here());
}

void VariableChange::print(std::ostream & os) const {}

}  // namespace ijedi

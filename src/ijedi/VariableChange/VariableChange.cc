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
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/State.h"
#include "mist/base/ModelData.h"
#include "oops/util/Logger.h"
#include "vader/vader.h"

namespace ijedi {

VariableChange::VariableChange(const eckit::Configuration &, const Geometry & geometry) {
  oops::Log::trace() << "ijedi::VariableChange::VariableChange starting" << std::endl;
  eckit::LocalConfiguration configCookbook{};
  const auto cb = ijedi::vaderCookbook();
  for (const auto & [key, val] : cb) {
    configCookbook.set(key, val);
  }

  auto configModelData = mist::base::ModelData(geometry).modelData();
  eckit::LocalConfiguration config{};
  config.set(vader::configCookbookKey, configCookbook);
  config.set(vader::configModelVarsKey, configModelData);

  varchange_ = std::make_unique<mist::utils::VariableChange>(config);
  oops::Log::trace() << "ijedi::VariableChange::VariableChange done" << std::endl;
}

void VariableChange::changeVar(State & xx, const oops::Variables & vars) const {
  oops::Log::trace() << "ijedi::VariableChange::changeVar starting" << std::endl;
  varchange_->changeVar(xx, vars);
  xx.setAtlasFieldMetadata();
  oops::Log::trace() << "ijedi::VariableChange::changeVar done" << std::endl;
}

void VariableChange::print(std::ostream & os) const {}

}  // namespace ijedi

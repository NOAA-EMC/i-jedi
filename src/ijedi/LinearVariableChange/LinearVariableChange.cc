/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/LinearVariableChange/LinearVariableChange.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/VariableChange/VaderCookbook.h"
#include "mist/base/ModelData.h"

namespace ijedi {

LinearVariableChange::LinearVariableChange(const Geometry & geometry,
                                           const eckit::Configuration & varchangeConfig) {
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

}  // namespace ijedi

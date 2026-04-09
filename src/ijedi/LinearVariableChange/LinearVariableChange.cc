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
                                           const eckit::Configuration &) {
  eckit::LocalConfiguration configCookbook{};
  const auto cb = ijedi::vaderCookbook();
  for (const auto & [key, val] : cb) {
    configCookbook.set(key, val);
  }

  auto configModelData = mist::base::ModelData(geometry).modelData();

  initVaderVariableChange(configCookbook, configModelData);
}

}  // namespace ijedi

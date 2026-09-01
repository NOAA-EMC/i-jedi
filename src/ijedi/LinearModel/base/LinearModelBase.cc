/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/LinearModel/base/LinearModelBase.h"

#include <map>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

LinearModelBase::LinearModelBase(const Geometry & /*geom*/, const eckit::Configuration & config)
  : tstep_(config.getString("time step")) {}

// -------------------------------------------------------------------------------------------------

LinearModelFactory::LinearModelFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in ijedi::LinearModelFactory."
                       << std::endl;
    ABORT("Element already registered in ijedi::LinearModelFactory.");
  }
  getMakers()[name] = this;
}

// -------------------------------------------------------------------------------------------------

LinearModelBase * LinearModelFactory::create(const Geometry & geom,
                                             const eckit::Configuration & config) {
  oops::Log::trace() << "ijedi::LinearModelFactory::create starting" << std::endl;

  const std::string id = config.getString("name");

  std::map<std::string, LinearModelFactory *>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    std::string registered;
    for (const std::string & name : getMakerNames()) {
      registered += (registered.empty() ? "" : ", ") + name;
    }
    throw eckit::BadValue("ijedi::LinearModelFactory: no linear model named '" + id +
                          "'. Registered linear models are: " + registered, Here());
  }

  LinearModelBase * ptr = jloc->second->make(geom, config);
  oops::Log::trace() << "ijedi::LinearModelFactory::create done" << std::endl;
  return ptr;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

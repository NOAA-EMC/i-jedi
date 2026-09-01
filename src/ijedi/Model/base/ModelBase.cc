/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/base/ModelBase.h"

#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

ModelBase::ModelBase(const Geometry & /*geom*/, const eckit::Configuration & config)
  : tstep_(config.getString("time step")), vars_() {
  if (config.has("model variables")) {
    vars_ = oops::Variables(config.getStringVector("model variables"));
  }
}

// -------------------------------------------------------------------------------------------------

ModelFactory::ModelFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in ijedi::ModelFactory." << std::endl;
    ABORT("Element already registered in ijedi::ModelFactory.");
  }
  getMakers()[name] = this;
}

// -------------------------------------------------------------------------------------------------

ModelBase * ModelFactory::create(const Geometry & geom, const eckit::Configuration & config) {
  oops::Log::trace() << "ijedi::ModelFactory::create starting" << std::endl;

  // "name" is the same key oops::ModelFactory dispatches on, so a single model block in the
  // yaml selects both the oops wrapper and the i-jedi implementation.
  const std::string id = config.getString("name");

  std::map<std::string, ModelFactory *>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    std::string registered;
    for (const std::string & name : getMakerNames()) {
      registered += (registered.empty() ? "" : ", ") + name;
    }
    throw eckit::BadValue("ijedi::ModelFactory: no model named '" + id +
                          "'. Registered models are: " + registered, Here());
  }

  ModelBase * ptr = jloc->second->make(geom, config);
  oops::Log::trace() << "ijedi::ModelFactory::create done" << std::endl;
  return ptr;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

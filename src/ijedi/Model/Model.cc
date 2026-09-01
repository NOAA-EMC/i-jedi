/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/Model.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

Model::Model(const Geometry & geom, const eckit::Configuration & config)
  : model_(ModelFactory::create(geom, config)) {
  oops::Log::trace() << "ijedi::Model created" << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

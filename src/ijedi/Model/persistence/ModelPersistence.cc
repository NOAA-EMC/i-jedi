/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/persistence/ModelPersistence.h"

#include <ostream>

#include "oops/util/Logger.h"

#include "ijedi/State/State.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

static ModelMaker<ModelPersistence> makerModelPersistence_("Persistence");

// -------------------------------------------------------------------------------------------------

void ModelPersistence::step(State & xx, const mist::ModelAuxControl & /*maux*/) const {
  oops::Log::trace() << "ModelPersistence::step starting" << std::endl;
  // The fields are left untouched; only the valid time moves. Advancing the time is not
  // optional: oops::Model::forecast loops until the state reaches the end of the leg.
  xx.updateTime(tstep_);
  oops::Log::trace() << "ModelPersistence::step done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void ModelPersistence::print(std::ostream & os) const {
  os << "ijedi::ModelPersistence, time step " << tstep_;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/LinearModel/persistence/LinearModelPersistence.h"

#include <ostream>

#include "oops/util/Logger.h"

#include "ijedi/Increment/Increment.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

static LinearModelMaker<LinearModelPersistence> makerLinearModelPersistence_("Persistence");

// -------------------------------------------------------------------------------------------------

void LinearModelPersistence::stepTL(Increment & dx, const mist::ModelAuxIncrement &) const {
  // Advancing the clock is not optional: oops::LinearModel::forecastTL loops until the
  // increment reaches the end of the window.
  dx.updateTime(tstep_);
}

// -------------------------------------------------------------------------------------------------

void LinearModelPersistence::stepAD(Increment & dx, mist::ModelAuxIncrement &) const {
  // The adjoint runs backwards through the window.
  dx.updateTime(-tstep_);
}

// -------------------------------------------------------------------------------------------------

void LinearModelPersistence::print(std::ostream & os) const {
  os << "ijedi::LinearModelPersistence, time step " << tstep_;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

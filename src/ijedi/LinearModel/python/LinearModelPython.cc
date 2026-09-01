/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/LinearModel/python/LinearModelPython.h"

#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/Model/python/FieldExchange.h"
#include "ijedi/State/State.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

static LinearModelMaker<LinearModelPython> makerLinearModelPython_("PythonModel");

// -------------------------------------------------------------------------------------------------

LinearModelPython::LinearModelPython(const Geometry & geom, const eckit::Configuration & config)
  : LinearModelBase(geom, config), geom_(geom),
    adapter_(config.getString("adapter")),
    backend_(config.getString("backend", "subprocess")) {
  oops::Log::trace() << "LinearModelPython::LinearModelPython starting" << std::endl;

  eckit::LocalConfiguration bridgeConf(config);
  bridgeConf.set("backend", backend_);
  bridgeConf.set("adapter", adapter_);
  bridge_ = PyModelBridge::create(bridgeConf, geom_.getComm());

  const PyModelBridge::Declaration & decl = bridge_->declaration();
  if (!decl.supportsLinear) {
    throw eckit::BadValue("LinearModelPython: the adapter '" + adapter_ + "' declares no "
                          "tangent linear or adjoint model. Use 'name: Identity' for the "
                          "linear model, or an adapter that provides one.", Here());
  }

  // Both the prognostic variables and the forcings, exactly as ModelPython does. The
  // trajectory is established by encoding a full state, and the adapter refuses to encode
  // without its forcings - so leaving them out of the map makes setTrajectory fail with a
  // complaint about a variable the linear model never appeared to ask for.
  for (const auto & entry : decl.inputVariables) {
    if (!entry.second.empty()) fieldNameMap_[entry.second] = entry.first;
  }
  for (const auto & entry : decl.forcingVariables) {
    if (!entry.second.empty()) fieldNameMap_[entry.second] = entry.first;
  }
  if (config.has("field name map")) {
    const eckit::LocalConfiguration mapConf(config, "field name map");
    for (const std::string & key : mapConf.keys()) {
      fieldNameMap_[key] = mapConf.getString(key);
    }
  }
  if (config.has("field scaling")) {
    const eckit::LocalConfiguration scaleConf(config, "field scaling");
    for (const std::string & key : scaleConf.keys()) {
      fieldScaling_[key] = scaleConf.getDouble(key);
    }
  }

  const util::Duration & inner = decl.timestep;
  if (inner.toSeconds() <= 0 || tstep_.toSeconds() % inner.toSeconds() != 0) {
    throw eckit::BadValue("LinearModelPython: 'time step' must be a positive whole multiple "
                          "of the adapter's internal timestep", Here());
  }
  innerSteps_ = static_cast<int>(tstep_.toSeconds() / inner.toSeconds());

  oops::Log::info() << "LinearModelPython: adapter '" << adapter_ << "', " << innerSteps_
                    << " internal steps per " << tstep_ << " linear step" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearModelPython::setTrajectory(const State & xx, State & /*xtraj*/,
                                      const mist::ModelAuxControl & /*maux*/) {
  oops::Log::trace() << "LinearModelPython::setTrajectory starting" << std::endl;

  PyModelBridge::Fields trajectory;
  gatherToGlobal(geom_.functionSpace(), xx.fieldSet(), fieldNameMap_, fieldScaling_, trajectory);
  if (!trajectory.empty()) {
    bridge_->setTrajectory(trajectory, xx.validTime());
  }

  oops::Log::trace() << "LinearModelPython::setTrajectory done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearModelPython::stepTL(Increment & dx, const mist::ModelAuxIncrement & /*maux*/) const {
  PyModelBridge::Fields in;
  gatherToGlobal(geom_.functionSpace(), dx.fieldSet(), fieldNameMap_, fieldScaling_, in);

  PyModelBridge::Fields out;
  if (geom_.getComm().rank() == 0) {
    bridge_->advanceTL(in, out, innerSteps_);
  }

  scatterFromGlobal(geom_.functionSpace(), out, fieldNameMap_, fieldScaling_, dx.fieldSet());
  dx.updateTime(tstep_);
}

// -------------------------------------------------------------------------------------------------

void LinearModelPython::stepAD(Increment & dx, mist::ModelAuxIncrement & /*maux*/) const {
  PyModelBridge::Fields in;
  gatherToGlobal(geom_.functionSpace(), dx.fieldSet(), fieldNameMap_, fieldScaling_, in);

  PyModelBridge::Fields out;
  if (geom_.getComm().rank() == 0) {
    bridge_->advanceAD(in, out, innerSteps_);
  }

  scatterFromGlobal(geom_.functionSpace(), out, fieldNameMap_, fieldScaling_, dx.fieldSet());
  // The adjoint runs backwards through the window.
  dx.updateTime(-tstep_);
}

// -------------------------------------------------------------------------------------------------

void LinearModelPython::print(std::ostream & os) const {
  os << "ijedi::LinearModelPython" << std::endl;
  os << "- adapter: " << adapter_ << " (" << bridge_->declaration().description << ")"
     << std::endl;
  os << "- linear step: " << tstep_ << " (" << innerSteps_ << " internal steps)";
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

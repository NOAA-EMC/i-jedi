/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/python/ModelPython.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Model/python/FieldExchange.h"
#include "ijedi/State/State.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

static ModelMaker<ModelPython> makerModelPython_("PythonModel");

// -------------------------------------------------------------------------------------------------

namespace {

// Pressure levels are compared in Pa; half a pascal separates any two ERA5 levels comfortably.
constexpr double kLevelTolerance = 0.5;

}  // namespace

// -------------------------------------------------------------------------------------------------

ModelPython::ModelPython(const Geometry & geom, const eckit::Configuration & config)
  : ModelBase(geom, config), geom_(geom),
    adapter_(config.getString("adapter")),
    backend_(config.getString("backend", "embedded")) {
  oops::Log::trace() << "ModelPython::ModelPython starting" << std::endl;

  // Bring the model up now rather than at the first forecast: a bad checkpoint or a geometry
  // that does not match the adapter should be a startup failure, not a surprise several
  // minutes into a run. Every task builds a bridge so that such a failure is raised
  // everywhere at once instead of deadlocking on one rank.
  eckit::LocalConfiguration bridgeConf(config);
  bridgeConf.set("backend", backend_);
  bridgeConf.set("adapter", adapter_);
  bridge_ = PyModelBridge::create(bridgeConf, geom_.getComm());

  const PyModelBridge::Declaration & decl = bridge_->declaration();
  checkDeclaration(decl);

  // The adapter's suggested JEDI names, overridable from the yaml. Inputs and forcings are
  // both exchanged, so they go into one map.
  for (const auto & entry : decl.inputVariables) {
    if (!entry.second.empty()) fieldNameMap_[entry.second] = entry.first;
  }
  for (const auto & entry : decl.forcingVariables) {
    if (!entry.second.empty()) fieldNameMap_[entry.second] = entry.first;
  }

  const std::vector<std::string> longNames = geom_.getFieldMetadata().getLongNames();
  const std::unordered_set<std::string> validNames(longNames.begin(), longNames.end());

  if (config.has("field name map")) {
    const eckit::LocalConfiguration mapConf(config, "field name map");
    for (const std::string & key : mapConf.keys()) {
      if (validNames.find(key) == validNames.end()) {
        throw eckit::BadValue("ModelPython: \"field name map\" contains \"" + key +
                              "\", which is not part of the field metadata.", Here());
      }
      fieldNameMap_[key] = mapConf.getString(key);
    }
  }

  if (config.has("field scaling")) {
    const eckit::LocalConfiguration scaleConf(config, "field scaling");
    for (const std::string & key : scaleConf.keys()) {
      if (validNames.find(key) == validNames.end()) {
        throw eckit::BadValue("ModelPython: \"field scaling\" contains \"" + key +
                              "\", which is not part of the field metadata.", Here());
      }
      fieldScaling_[key] = scaleConf.getDouble(key);
    }
  }

  // Every JEDI name we intend to exchange must be a name the metadata knows, or the state
  // cannot carry it.
  for (const auto & entry : fieldNameMap_) {
    if (validNames.find(entry.first) == validNames.end()) {
      throw eckit::BadValue("ModelPython: the adapter maps its variable '" + entry.second +
                            "' onto JEDI long name '" + entry.first + "', which is not part "
                            "of the field metadata. Override it with \"field name map\".",
                            Here());
    }
  }

  // Every variable the adapter declared has to be reachable from some JEDI name.
  std::unordered_set<std::string> mappedAdapterNames;
  for (const auto & entry : fieldNameMap_) mappedAdapterNames.insert(entry.second);
  for (const auto & entry : decl.inputVariables) {
    if (mappedAdapterNames.find(entry.first) == mappedAdapterNames.end()) {
      throw eckit::BadValue("ModelPython: the adapter needs '" + entry.first + "' but nothing "
                            "maps onto it. Add it to \"field name map\".", Here());
    }
  }

  // How many internal steps make up one oops step.
  const util::Duration & inner = decl.timestep;
  if (inner.toSeconds() <= 0) {
    throw eckit::BadValue("ModelPython: the adapter reports a non-positive timestep", Here());
  }
  if (tstep_.toSeconds() % inner.toSeconds() != 0) {
    std::stringstream msg;
    msg << "ModelPython: 'time step' " << tstep_ << " is not a whole multiple of the adapter's "
        << "internal timestep " << inner;
    throw eckit::BadValue(msg.str(), Here());
  }
  innerSteps_ = static_cast<int>(tstep_.toSeconds() / inner.toSeconds());

  // Unless the yaml says otherwise, the model variables are the ones it exchanges.
  if (vars_.size() == 0) {
    std::vector<std::string> names;
    for (const auto & entry : fieldNameMap_) names.push_back(entry.first);
    vars_ = oops::Variables(names);
  }

  oops::Log::info() << "ModelPython: adapter '" << adapter_ << "' (" << decl.description
                    << ") via backend '" << backend_ << "', " << innerSteps_
                    << " internal steps per " << tstep_ << " output step" << std::endl;

  oops::Log::trace() << "ModelPython::ModelPython done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void ModelPython::checkDeclaration(const PyModelBridge::Declaration & decl) {
  if (decl.levelsPa.empty()) return;

  if (!geom_.modelData().has("vertical_coordinate_reference_pressure")) {
    throw eckit::BadValue("ModelPython: the adapter is discretised on pressure levels but the "
                          "geometry carries no reference pressure column. Set 'vertical "
                          "levels' on the geometry to the adapter's levels.", Here());
  }

  const std::vector<double> geomLevels =
      geom_.modelData().getDoubleVector("vertical_coordinate_reference_pressure");

  if (geomLevels.size() != decl.levelsPa.size()) {
    std::stringstream msg;
    msg << "ModelPython: the adapter declares " << decl.levelsPa.size() << " pressure levels "
        << "but the geometry has " << geomLevels.size();
    throw eckit::BadValue(msg.str(), Here());
  }

  // Compare values, not just counts: two 37-level configurations that disagree about which
  // 37 levels would otherwise run and produce nonsense.
  for (size_t k = 0; k < geomLevels.size(); ++k) {
    if (std::abs(geomLevels[k] - decl.levelsPa[k]) > kLevelTolerance) {
      std::stringstream msg;
      msg << "ModelPython: at level " << k << " the adapter expects " << decl.levelsPa[k]
          << " Pa but the geometry has " << geomLevels[k] << " Pa";
      throw eckit::BadValue(msg.str(), Here());
    }
  }
}

// -------------------------------------------------------------------------------------------------

void ModelPython::initialize(State & xx) const {
  oops::Log::trace() << "ModelPython::initialize starting" << std::endl;

  PyModelBridge::Fields inputs;
  gatherToGlobal(geom_.functionSpace(), xx.fieldSet(), fieldNameMap_, fieldScaling_, inputs);

  // Populated on the root task only, which is the task whose bridge does the work.
  if (!inputs.empty()) {
    bridge_->encode(inputs, xx.validTime());
  }

  oops::Log::trace() << "ModelPython::initialize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void ModelPython::step(State & xx, const mist::ModelAuxControl & /*maux*/) const {
  oops::Log::trace() << "ModelPython::step starting" << std::endl;

  PyModelBridge::Fields outputs;
  if (geom_.getComm().rank() == 0) {
    bridge_->advance(innerSteps_);
    bridge_->decode(outputs);
  }

  scatterFromGlobal(geom_.functionSpace(), outputs, fieldNameMap_, fieldScaling_,
                    xx.fieldSet());

  xx.updateTime(tstep_);

  oops::Log::trace() << "ModelPython::step done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void ModelPython::finalize(State & /*xx*/) const {
  oops::Log::trace() << "ModelPython::finalize starting" << std::endl;
  bridge_->reset();
  oops::Log::trace() << "ModelPython::finalize done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void ModelPython::print(std::ostream & os) const {
  os << "ijedi::ModelPython" << std::endl;
  os << "- adapter: " << adapter_ << " (" << bridge_->declaration().description << ")"
     << std::endl;
  os << "- backend: " << backend_ << std::endl;
  os << "- output interval: " << tstep_ << " (" << innerSteps_ << " internal steps)"
     << std::endl;
  os << "- linear model: " << (bridge_->declaration().supportsLinear ? "yes" : "no")
     << std::endl;
  os << "- variables exchanged:";
  for (const auto & entry : fieldNameMap_) {
    os << std::endl << "    " << entry.first << " -> " << entry.second;
  }
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

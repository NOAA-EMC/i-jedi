/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ijedi/LinearModel/LinearModelBase.h"

namespace mist {
class ModelAuxControl;
class ModelAuxIncrement;
}  // namespace mist

namespace ijedi {

class Geometry;
class Increment;
class State;

// Dispatches to the linear model registered in ijedi::LinearModelFactory under the YAML
// `name:` (e.g. "SQG"). oops registers one linear model per entry of names(); the generic
// linear models oops provides itself (e.g. "Identity") remain available as before.
class LinearModel : public util::Printable, private util::ObjectCounter<LinearModel> {
 public:
  static const std::string classname() { return "ijedi::LinearModel"; }
  static std::vector<std::string> names() { return LinearModelFactory::getMakerNames(); }

  LinearModel(const Geometry & geom, const eckit::Configuration & config):
    tstep_(config.getString("time step")),
    impl_(LinearModelFactory::create(config.getString("name"), geom, config, tstep_)) {}
  ~LinearModel() = default;

  /// Model trajectory computation
  void setTrajectory(const State & xx, State &, const mist::ModelAuxControl &) {
    impl_->setTrajectory(xx);
  }

  /// Run TLM and its adjoint (stepTL/stepAD advance/rewind the increment's valid time)
  void initializeTL(Increment & dx) const { impl_->initializeTL(dx); }
  void stepTL(Increment & dx, const mist::ModelAuxIncrement &) const { impl_->stepTL(dx); }
  void finalizeTL(Increment & dx) const { impl_->finalizeTL(dx); }

  void initializeAD(Increment & dx) const { impl_->initializeAD(dx); }
  void stepAD(Increment & dx, mist::ModelAuxIncrement &) const { impl_->stepAD(dx); }
  void finalizeAD(Increment & dx) const { impl_->finalizeAD(dx); }

  /// Other utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const util::Duration & stepTrajectory() const {return tstep_;}

 private:
  void print(std::ostream & os) const override { impl_->print(os); }

  util::Duration tstep_;
  std::unique_ptr<LinearModelBase> impl_;
};

}  // namespace ijedi

/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace mist::base {
class ModelAuxControl;
class ModelAuxIncrement;
}  // namespace mist::base

namespace ijedi {

class Geometry;
class Increment;
class State;

// An identity linear model.
class LinearModel : public util::Printable, private util::ObjectCounter<LinearModel> {
 public:
  static const std::string classname() { return "ijedi::Model"; }
  static std::vector<std::string> names() { return {}; }

  LinearModel(const Geometry &, const eckit::Configuration & config):
    tstep_(config.getString("time step")) {}
  ~LinearModel() = default;

  /// Model trajectory computation
  void setTrajectory(const State &, State &, const mist::base::ModelAuxControl &) {}

/// Run TLM and its adjoint
  void initializeTL(Increment &) const {}
  void stepTL(Increment &, const mist::base::ModelAuxIncrement &) const {}
  void finalizeTL(Increment &) const {}

  void initializeAD(Increment &) const {}
  void stepAD(Increment &, mist::base::ModelAuxIncrement &) const {}
  void finalizeAD(Increment &) const {}

/// Other utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const util::Duration & stepTrajectory() const {return tstep_;}

 private:
  void print(std::ostream &) const override {}

  util::Duration tstep_;
};

}  // namespace ijedi

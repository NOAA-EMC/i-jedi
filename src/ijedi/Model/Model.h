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

namespace mist {
class ModelAuxControl;
}  // namespace mist

namespace ijedi {

class Geometry;
class State;

// A persistence model.
class Model : public util::Printable, private util::ObjectCounter<Model> {
 public:
  static const std::string classname() { return "ijedi::Model"; }
  static std::vector<std::string> names() { return {}; }

  Model(const Geometry &, const eckit::Configuration & config):
    tstep_(config.getString("time step")), vars_() {}
  ~Model() = default;

  void initialize(State &) const {}
  void step(State &, const mist::ModelAuxControl &) const {}
  void finalize(State &) const {}

  const util::Duration & timeResolution() const { return tstep_; }
  const oops::Variables & variables() const { return vars_; }

 private:
  void print(std::ostream &) const override {}

  util::Duration tstep_;
  const oops::Variables vars_;
};

}  // namespace ijedi

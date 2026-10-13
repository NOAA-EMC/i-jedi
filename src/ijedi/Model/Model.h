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

#include "ijedi/Model/ModelBase.h"
#include "ijedi/State/State.h"

namespace mist {
class ModelAuxControl;
}  // namespace mist

namespace ijedi {

class Geometry;

// Dispatches to the forecast model registered in ijedi::ModelFactory under the YAML
// `name:` (e.g. "SQG"). oops registers one model per entry of names(); the generic models
// oops provides itself ("Identity", "PseudoModel") remain available as before.
class Model : public util::Printable, private util::ObjectCounter<Model> {
 public:
  static const std::string classname() { return "ijedi::Model"; }
  static std::vector<std::string> names() { return ModelFactory::getMakerNames(); }

  Model(const Geometry & geom, const eckit::Configuration & config):
    tstep_(config.getString("time step")),
    impl_(ModelFactory::create(config.getString("name"), geom, config, tstep_)) {}
  ~Model() = default;

  void initialize(State & xx) const { impl_->initialize(xx); }
  // oops' forecast loop relies on step() advancing the state's valid time.
  void step(State & xx, const mist::ModelAuxControl &) const {
    impl_->step(xx);
    xx.updateTime(tstep_);
  }
  void finalize(State & xx) const { impl_->finalize(xx); }

  const util::Duration & timeResolution() const { return tstep_; }
  const oops::Variables & variables() const { return impl_->variables(); }

 private:
  void print(std::ostream & os) const override { impl_->print(os); }

  util::Duration tstep_;
  std::unique_ptr<ModelBase> impl_;
};

}  // namespace ijedi

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

#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ijedi/Model/base/ModelBase.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace mist {
class ModelAuxControl;
}  // namespace mist

namespace ijedi {

class Geometry;
class State;

// -------------------------------------------------------------------------------------------------

/// \brief The MODEL-specific forecast model oops sees, dispatching to one of i-jedi's model
///        implementations.
///
/// oops::instantiateModelFactory registers one oops::interface::Model per entry in names(),
/// and oops::ModelFactory then selects between them on the "name" key. names() returns the
/// implementations registered with ijedi::ModelFactory, so the two factories agree by
/// construction and a newly added model becomes reachable from yaml without touching this
/// class.
class Model : public util::Printable, private util::ObjectCounter<Model> {
 public:
  static const std::string classname() { return "ijedi::Model"; }
  static std::vector<std::string> names() { return ModelFactory::getMakerNames(); }

  Model(const Geometry &, const eckit::Configuration &);
  ~Model() = default;

  void initialize(State & xx) const { model_->initialize(xx); }
  void step(State & xx, const mist::ModelAuxControl & maux) const { model_->step(xx, maux); }
  void finalize(State & xx) const { model_->finalize(xx); }

  const util::Duration & timeResolution() const { return model_->timeResolution(); }
  const oops::Variables & variables() const { return model_->variables(); }

 private:
  void print(std::ostream & os) const override { os << *model_; }

  std::unique_ptr<ModelBase> model_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

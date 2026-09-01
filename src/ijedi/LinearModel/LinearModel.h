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

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ijedi/LinearModel/base/LinearModelBase.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace mist {
class ModelAuxControl;
class ModelAuxIncrement;
}  // namespace mist

namespace ijedi {

class Geometry;
class Increment;
class State;

// -------------------------------------------------------------------------------------------------

/// \brief The MODEL-specific linear model oops sees, dispatching to one of i-jedi's
///        implementations. The same arrangement as ijedi::Model, and for the same reasons.
class LinearModel : public util::Printable, private util::ObjectCounter<LinearModel> {
 public:
  static const std::string classname() { return "ijedi::LinearModel"; }
  static std::vector<std::string> names() { return LinearModelFactory::getMakerNames(); }

  LinearModel(const Geometry &, const eckit::Configuration &);
  ~LinearModel() = default;

  void setTrajectory(const State & xx, State & xtraj, const mist::ModelAuxControl & maux)
    { model_->setTrajectory(xx, xtraj, maux); }

  void initializeTL(Increment & dx) const { model_->initializeTL(dx); }
  void stepTL(Increment & dx, const mist::ModelAuxIncrement & maux) const
    { model_->stepTL(dx, maux); }
  void finalizeTL(Increment & dx) const { model_->finalizeTL(dx); }

  void initializeAD(Increment & dx) const { model_->initializeAD(dx); }
  void stepAD(Increment & dx, mist::ModelAuxIncrement & maux) const
    { model_->stepAD(dx, maux); }
  void finalizeAD(Increment & dx) const { model_->finalizeAD(dx); }

  const util::Duration & timeResolution() const { return model_->timeResolution(); }
  const util::Duration & stepTrajectory() const { return model_->stepTrajectory(); }

 private:
  void print(std::ostream & os) const override { os << *model_; }

  std::unique_ptr<LinearModelBase> model_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "ijedi/LinearModel/base/LinearModelBase.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace ijedi {

class Geometry;

// -------------------------------------------------------------------------------------------------

/// \brief Persistence: perturbations are carried unchanged, only the clock moves.
///
/// This is its own adjoint, so it passes the dot-product test trivially and is exactly linear.
/// That makes it a clean check of everything around it - i-jedi's linear model factory, the
/// dispatcher, the increment handling and the forward and backward clock - and useless as a
/// statement about any model's actual dynamics.
///
/// Registered as "Persistence" rather than "Identity": oops::instantiateLinearModelFactory
/// registers its own generic IdentityLinearModel under that name, and two makers cannot share
/// one. The name also matches ijedi's nonlinear ModelPersistence, which does the same thing.
class LinearModelPersistence : public LinearModelBase {
 public:
  static const std::string classname() { return "ijedi::LinearModelPersistence"; }

  LinearModelPersistence(const Geometry & geom, const eckit::Configuration & config)
    : LinearModelBase(geom, config) {}
  ~LinearModelPersistence() = default;

  void setTrajectory(const State &, State &, const mist::ModelAuxControl &) override {}

  void initializeTL(Increment &) const override {}
  void stepTL(Increment &, const mist::ModelAuxIncrement &) const override;
  void finalizeTL(Increment &) const override {}

  void initializeAD(Increment &) const override {}
  void stepAD(Increment &, mist::ModelAuxIncrement &) const override;
  void finalizeAD(Increment &) const override {}

 private:
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <ostream>
#include <string>

#include "ijedi/LinearModel/base/LinearModelBase.h"
#include "ijedi/Model/python/PyModelBridge.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace ijedi {

class Geometry;
class Increment;
class State;

// -------------------------------------------------------------------------------------------------

/// \brief The tangent linear and adjoint of a Python forecast model.
///
/// This is the reason a differentiable Python model is interesting for data assimilation.
/// Both JAX and PyTorch derive forward- and reverse-mode derivatives from the forward code, so
/// what is a multi-year hand-coding exercise for a conventional model is a wrapping exercise
/// here - the adapter calls jax.jvp and jax.vjp, and the result is an adjoint that is the
/// exact transpose of the tangent linear rather than an approximation maintained by hand.
///
/// Two properties of the linearisation are the adapter's responsibility, and both are traps
/// that the usual dot-product test cannot detect, because jvp and vjp are exact transposes of
/// whatever function they are handed: the operator must be the whole composition
/// decode(advance(encode(x))), since the encoder and decoder are themselves learned; and for a
/// stochastic model the random draw must be the one the trajectory used.
class LinearModelPython : public LinearModelBase {
 public:
  static const std::string classname() { return "ijedi::LinearModelPython"; }

  LinearModelPython(const Geometry &, const eckit::Configuration &);
  ~LinearModelPython() = default;

  void setTrajectory(const State &, State &, const mist::ModelAuxControl &) override;

  void initializeTL(Increment &) const override {}
  void stepTL(Increment &, const mist::ModelAuxIncrement &) const override;
  void finalizeTL(Increment &) const override {}

  void initializeAD(Increment &) const override {}
  void stepAD(Increment &, mist::ModelAuxIncrement &) const override;
  void finalizeAD(Increment &) const override {}

 private:
  void print(std::ostream &) const override;

  const Geometry & geom_;
  std::string adapter_;
  std::string backend_;
  std::map<std::string, std::string> fieldNameMap_;
  std::map<std::string, double> fieldScaling_;

  std::unique_ptr<PyModelBridge> bridge_;
  int innerSteps_ = 0;
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

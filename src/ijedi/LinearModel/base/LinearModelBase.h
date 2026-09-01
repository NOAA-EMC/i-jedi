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
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/util/AssociativeContainers.h"
#include "oops/util/Duration.h"
#include "oops/util/Printable.h"

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

/// \brief Base class for i-jedi's tangent linear and adjoint models.
///
/// The same polymorphism as ModelBase, and for the same reason: i-jedi supports several
/// models at once, so it has to support several linear models. The concrete one is chosen at
/// runtime from the "name" key, which is what oops::LinearModelFactory dispatches on.
///
/// A linear model derived automatically from a differentiable forecast model is the reason
/// this is worth building at all. For a conventional model, writing a tangent linear and its
/// adjoint is a multi-year hand-coding effort; for a model written in JAX or PyTorch it is a
/// wrapping exercise, because the framework differentiates the forward code for you.
class LinearModelBase : public util::Printable, private boost::noncopyable {
 public:
  LinearModelBase(const Geometry &, const eckit::Configuration &);
  virtual ~LinearModelBase() = default;

  /// \brief Establish the nonlinear trajectory the linear operators linearise about. Called
  ///        once per window before any TL or AD step.
  virtual void setTrajectory(const State &, State &, const mist::ModelAuxControl &) = 0;

  virtual void initializeTL(Increment &) const = 0;
  virtual void stepTL(Increment &, const mist::ModelAuxIncrement &) const = 0;
  virtual void finalizeTL(Increment &) const = 0;

  virtual void initializeAD(Increment &) const = 0;
  virtual void stepAD(Increment &, mist::ModelAuxIncrement &) const = 0;
  virtual void finalizeAD(Increment &) const = 0;

  const util::Duration & timeResolution() const { return tstep_; }
  const util::Duration & stepTrajectory() const { return tstep_; }

 protected:
  util::Duration tstep_;
};

// -------------------------------------------------------------------------------------------------

class LinearModelFactory {
 public:
  static LinearModelBase * create(const Geometry &, const eckit::Configuration &);

  static std::vector<std::string> getMakerNames() { return oops::keys(getMakers()); }

  virtual ~LinearModelFactory() = default;

 protected:
  explicit LinearModelFactory(const std::string & name);

 private:
  virtual LinearModelBase * make(const Geometry &, const eckit::Configuration &) = 0;

  static std::map<std::string, LinearModelFactory *> & getMakers() {
    static std::map<std::string, LinearModelFactory *> makers_;
    return makers_;
  }
};

// -------------------------------------------------------------------------------------------------

template <class T>
class LinearModelMaker : public LinearModelFactory {
  LinearModelBase * make(const Geometry & geom, const eckit::Configuration & config) override {
    return new T(geom, config);
  }

 public:
  explicit LinearModelMaker(const std::string & name) : LinearModelFactory(name) {}
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

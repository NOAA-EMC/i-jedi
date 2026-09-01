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

#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Duration.h"
#include "oops/util/Printable.h"

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

/// \brief Base class for the forecast models i-jedi can drive.
///
/// i-jedi supports several native grids concurrently, so it also has to support several
/// forecast models. This mirrors the polymorphism already used for Geometry and Io: the
/// concrete model is chosen at runtime from the "name" key of the model configuration, which
/// is the same key oops::ModelFactory dispatches on.
///
/// The three-phase lifecycle is oops's, not i-jedi's: initialize() is called once before a
/// forecast leg, step() once per timeResolution() until the leg is over, and finalize() once
/// at the end. Models that carry expensive internal state should build it in initialize()
/// and release it in finalize() rather than rebuilding it every step.
class ModelBase : public util::Printable, private boost::noncopyable {
 public:
  ModelBase(const Geometry &, const eckit::Configuration &);
  virtual ~ModelBase() = default;

  virtual void initialize(State &) const = 0;
  virtual void step(State &, const mist::ModelAuxControl &) const = 0;
  virtual void finalize(State &) const = 0;

  /// \brief Interval at which oops sees the state, i.e. the model's output frequency. For a
  ///        model with its own shorter internal timestep this is not that timestep.
  const util::Duration & timeResolution() const { return tstep_; }

  /// \brief Variables the model needs to be handed, and hands back.
  const oops::Variables & variables() const { return vars_; }

 protected:
  util::Duration tstep_;
  oops::Variables vars_;
};

// -------------------------------------------------------------------------------------------------

class ModelFactory {
 public:
  static ModelBase * create(const Geometry &, const eckit::Configuration &);

  static std::vector<std::string> getMakerNames() { return oops::keys(getMakers()); }

  virtual ~ModelFactory() = default;

 protected:
  explicit ModelFactory(const std::string & name);

 private:
  virtual ModelBase * make(const Geometry &, const eckit::Configuration &) = 0;

  static std::map<std::string, ModelFactory *> & getMakers() {
    static std::map<std::string, ModelFactory *> makers_;
    return makers_;
  }
};

// -------------------------------------------------------------------------------------------------

template <class T>
class ModelMaker : public ModelFactory {
  ModelBase * make(const Geometry & geom, const eckit::Configuration & config) override {
    return new T(geom, config);
  }

 public:
  explicit ModelMaker(const std::string & name) : ModelFactory(name) {}
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

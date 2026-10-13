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

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/AssociativeContainers.h"
#include "oops/util/Duration.h"

namespace ijedi {

class Geometry;
class Increment;
class State;

// -------------------------------------------------------------------------------------------------
/// Base class for the linear (tangent linear and adjoint) models that ijedi::LinearModel
/// dispatches to, registered with LinearModelFactory under the YAML `linear model: name:`.
/// stepTL must advance the increment's valid time by one time step and stepAD move it back.
class LinearModelBase {
 public:
  virtual ~LinearModelBase() = default;

  virtual void setTrajectory(const State &) = 0;

  virtual void initializeTL(Increment &) const {}
  virtual void stepTL(Increment &) const = 0;
  virtual void finalizeTL(Increment &) const {}

  virtual void initializeAD(Increment &) const {}
  virtual void stepAD(Increment &) const = 0;
  virtual void finalizeAD(Increment &) const {}

  virtual void print(std::ostream &) const = 0;
};

// -------------------------------------------------------------------------------------------------
class LinearModelFactory {
 public:
  static std::unique_ptr<LinearModelBase> create(const std::string & name, const Geometry &,
                                                 const eckit::Configuration &,
                                                 const util::Duration & tstep);
  static std::vector<std::string> getMakerNames() { return oops::keys(getMakers()); }
  virtual ~LinearModelFactory() = default;

 protected:
  explicit LinearModelFactory(const std::string & name) {
    if (getMakers().find(name) != getMakers().end()) {
      throw eckit::BadValue(name + " already registered in ijedi::LinearModelFactory.", Here());
    }
    getMakers()[name] = this;
  }

 private:
  virtual LinearModelBase * make(const Geometry &, const eckit::Configuration &,
                                 const util::Duration &) = 0;
  static std::map<std::string, LinearModelFactory *> & getMakers() {
    static std::map<std::string, LinearModelFactory *> makers_;
    return makers_;
  }
};

// -------------------------------------------------------------------------------------------------
template <class T>
class LinearModelMaker : public LinearModelFactory {
  LinearModelBase * make(const Geometry & geom, const eckit::Configuration & conf,
                         const util::Duration & tstep) override {
    return new T(geom, conf, tstep);
  }

 public:
  explicit LinearModelMaker(const std::string & name) : LinearModelFactory(name) {}
};

// -------------------------------------------------------------------------------------------------
inline std::unique_ptr<LinearModelBase> LinearModelFactory::create(
    const std::string & name, const Geometry & geom, const eckit::Configuration & conf,
    const util::Duration & tstep) {
  auto it = getMakers().find(name);
  if (it == getMakers().end()) {
    throw eckit::BadValue("ijedi::LinearModelFactory: no linear model named \"" + name + "\"",
                          Here());
  }
  return std::unique_ptr<LinearModelBase>(it->second->make(geom, conf, tstep));
}

}  // namespace ijedi

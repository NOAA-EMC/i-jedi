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

#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Duration.h"

namespace ijedi {

class Geometry;
class State;

// -------------------------------------------------------------------------------------------------
/// Base class for the forecast models that ijedi::Model dispatches to. Each implementation
/// registers itself with ModelFactory under the name used as `model: name:` in the YAML.
class ModelBase {
 public:
  virtual ~ModelBase() = default;

  virtual void initialize(State &) const {}
  virtual void step(State &) const = 0;
  virtual void finalize(State &) const {}

  virtual const oops::Variables & variables() const = 0;
  virtual void print(std::ostream &) const = 0;
};

// -------------------------------------------------------------------------------------------------
class ModelFactory {
 public:
  static std::unique_ptr<ModelBase> create(const std::string & name, const Geometry &,
                                           const eckit::Configuration &,
                                           const util::Duration & tstep);
  static std::vector<std::string> getMakerNames() { return oops::keys(getMakers()); }
  virtual ~ModelFactory() = default;

 protected:
  explicit ModelFactory(const std::string & name) {
    if (getMakers().find(name) != getMakers().end()) {
      throw eckit::BadValue(name + " already registered in ijedi::ModelFactory.", Here());
    }
    getMakers()[name] = this;
  }

 private:
  virtual ModelBase * make(const Geometry &, const eckit::Configuration &,
                           const util::Duration &) = 0;
  static std::map<std::string, ModelFactory *> & getMakers() {
    static std::map<std::string, ModelFactory *> makers_;
    return makers_;
  }
};

// -------------------------------------------------------------------------------------------------
template <class T>
class ModelMaker : public ModelFactory {
  ModelBase * make(const Geometry & geom, const eckit::Configuration & conf,
                   const util::Duration & tstep) override {
    return new T(geom, conf, tstep);
  }

 public:
  explicit ModelMaker(const std::string & name) : ModelFactory(name) {}
};

// -------------------------------------------------------------------------------------------------
inline std::unique_ptr<ModelBase> ModelFactory::create(const std::string & name,
                                                       const Geometry & geom,
                                                       const eckit::Configuration & conf,
                                                       const util::Duration & tstep) {
  auto it = getMakers().find(name);
  if (it == getMakers().end()) {
    throw eckit::BadValue("ijedi::ModelFactory: no model named \"" + name + "\"", Here());
  }
  return std::unique_ptr<ModelBase>(it->second->make(geom, conf, tstep));
}

}  // namespace ijedi

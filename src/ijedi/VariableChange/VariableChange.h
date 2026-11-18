/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>

#include "mist/base/VariableChange.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace oops {
class Variables;
}  // namespace oops

namespace ijedi {

class Geometry;
class State;

class VariableChange : public mist::VariableChange {
 public:
  VariableChange(const eckit::Configuration &, const Geometry &);
  ~VariableChange() = default;

  // Perform transforms
  void changeVar(State &, const oops::Variables &) const;
  void changeVarInverse(State &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;
  const Geometry & geom_;
};

}  // namespace ijedi

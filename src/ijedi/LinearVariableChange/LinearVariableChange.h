/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "mist/base/LinearVariableChange.h"

namespace oops {
class Variables;
}  // namespace oops

namespace ijedi {

class Geometry;
class State;

class LinearVariableChange : public mist::base::LinearVariableChange {
 public:
  LinearVariableChange(const Geometry &, const eckit::Configuration &);

  /// Inject geometry-sourced trajectory ingredients (latitude, longitude,
  /// sea_area_fraction) required by SeaWaterTemperature_B's Jacobian before
  /// handing the trajectory to the base class.
  void changeVarTraj(const State &, const oops::Variables &);

 private:
  const Geometry & geom_;
};

}  // namespace ijedi

// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <string>

#include "ijedi/ErrorCovariance/ErrorCovariance.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/LinearVariableChange/LinearVariableChange.h"
#include "ijedi/ModelData.h"
#include "ijedi/State/State.h"

namespace ijedi {

struct Traits {
  static std::string name() {return "ijedi";}
  static std::string nameCovar() {return "ijediError";}

  typedef ijedi::ErrorCovariance       Covariance;
  typedef ijedi::Geometry              Geometry;
  typedef ijedi::Increment             Increment;
  typedef ijedi::LinearVariableChange  LinearVariableChange;
  typedef ijedi::ModelData             ModelData;
  typedef ijedi::State                 State;
};

}  // namespace ijedi

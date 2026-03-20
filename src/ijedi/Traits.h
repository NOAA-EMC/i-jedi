// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <string>

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/State.h"

namespace ijedi {

struct Traits {
  static std::string name() {return "ijedi";}
  static std::string nameCovar() {return "ijediError";}

  typedef ijedi::Geometry           Geometry;
  typedef ijedi::State              State;
};

}  // namespace ijedi

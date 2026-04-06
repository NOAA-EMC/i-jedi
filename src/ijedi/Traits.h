#pragma once

#include <string>

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/State.h"

namespace ijedi
{

  struct Traits
  {
    static std::string name() { return "ijedi"; }
    static std::string nameCovar() { return "ijediError"; }

    typedef ijedi::Geometry Geometry;
    typedef ijedi::State    State;
  };

}  // namespace ijedi

#pragma once

#include <string>

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/State/State.h"

#include "mist/base/ModelData.h"
#include "mist/base/VariableChange.h"

#include "oops/generic/UnstructuredInterpolator.h"

namespace ijedi
{

  struct Traits
  {
    static std::string name() { return "ijedi"; }
    static std::string nameCovar() { return "ijediError"; }

    typedef ijedi::Geometry                Geometry;
    typedef ijedi::State                   State;
    typedef mist::base::ModelData          ModelData;
    typedef mist::base::VariableChange     VariableChange;
    typedef ijedi::Increment               Increment;
    typedef oops::UnstructuredInterpolator LocalInterpolator;
  };

}  // namespace ijedi

// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace ijedi
{

  class Geometry;

  // ---------------------------------------------------------------------------
  // Stub ModelData required by SABER's SaberOuterBlockChain.
  // Returns empty model data — ijedi has no model-specific data to expose.
  // ---------------------------------------------------------------------------
  class ModelData : public util::Printable
  {
   public:
    static const std::string classname() { return "ijedi::ModelData"; }

    explicit ModelData(const Geometry & /*geom*/) {}
    ~ModelData() = default;

    static oops::Variables defaultVariables() { return oops::Variables(); }

    eckit::LocalConfiguration modelData() const
    {
      return eckit::LocalConfiguration();
    }

   private:
    void print(std::ostream & os) const override
    {
      os << "ModelData: empty stub";
    }
  };

}  // namespace ijedi

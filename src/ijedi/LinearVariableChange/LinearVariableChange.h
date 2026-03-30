// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class Geometry;
  class Increment;
  class State;

  // ---------------------------------------------------------------------------
  // Identity (no-op) linear variable change.
  //
  // Required by oops::ModelSpaceCovarianceBase<MODEL> even when no variable
  // change is needed.  All five changeVar* methods are no-ops.
  // ---------------------------------------------------------------------------
  class LinearVariableChange : public util::Printable
  {
   public:
    static const std::string classname() { return "ijedi::LinearVariableChange"; }

    LinearVariableChange(const Geometry &, const eckit::Configuration &);
    ~LinearVariableChange() = default;

    void changeVarTraj(const State &, const oops::Variables &) {}
    void changeVarTL(Increment &, const oops::Variables &) const {}
    void changeVarInverseTL(Increment &, const oops::Variables &) const {}
    void changeVarAD(Increment &, const oops::Variables &) const {}
    void changeVarInverseAD(Increment &, const oops::Variables &) const {}

   private:
    void print(std::ostream & os) const override
    {
      os << "LinearVariableChange: identity";
    }
  };

}  // namespace ijedi

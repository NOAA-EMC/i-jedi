// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <ostream>

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/base/StateBase.h"

namespace ijedi
{

  // ---------------------------------------------------------------------------
  // MOM6 concrete State implementation.
  //
  // Constructors call the appropriate StateBase::init* helper to set
  // functionSpace_, numLevels_, vars_, time_, and allocate fields.
  // All field operations (zero, accumul, norm, serialize, print, ...) are
  // inherited from StateBase.  Only read and write are model-specific.
  // ---------------------------------------------------------------------------
  class StateMOM6 : public StateBase
  {
   public:
    StateMOM6(const Geometry &, const oops::Variables &,
              const util::DateTime &);
    StateMOM6(const Geometry &, const eckit::Configuration &);
    StateMOM6(const Geometry &, const StateBase &);     // change geometry
    StateMOM6(const oops::Variables &, const StateMOM6 &);  // change vars

    void read(const eckit::Configuration &) override;
    void write(const eckit::Configuration &) const override;

   private:
    void fillAnalytic_();
    const Geometry & geom_;
  };

}  // namespace ijedi

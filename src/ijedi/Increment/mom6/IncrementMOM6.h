// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <ostream>

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Increment/base/IncrementBase.h"

namespace ijedi
{

  // ---------------------------------------------------------------------------
  // MOM6 concrete Increment implementation.
  //
  // Constructors call the appropriate IncrementBase::init* helper to set
  // functionSpace_, numLevels_, vars_, time_, and allocate fields.
  // All field operations are inherited from IncrementBase.
  // Only read and write are model-specific.
  // ---------------------------------------------------------------------------
  class IncrementMOM6 : public IncrementBase
  {
   public:
    IncrementMOM6(const Geometry &, const oops::Variables &,
                  const util::DateTime &);
    IncrementMOM6(const Geometry &, const IncrementBase &);  // change geom
    IncrementMOM6(const IncrementMOM6 &, bool copy = true);  // copy/nocopy

    void read(const eckit::Configuration &) override;
    void write(const eckit::Configuration &) const override;

   private:
    const Geometry & geom_;
  };

}  // namespace ijedi

// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <ostream>

#include "eckit/mpi/Comm.h"

#include "atlas/grid.h"

#include "ijedi/Geometry/base/GeometryBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class GeometryMPAS : public GeometryBase
  {
   public:
    GeometryMPAS(const eckit::Configuration &, const eckit::mpi::Comm &);
    void print(std::ostream &) const override;

    // Unified access to grid-specific parameters
    eckit::LocalConfiguration gridSpecific() const override;
  };

}  // namespace ijedi

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

  class GeometryFV3 : public GeometryBase
  {
   public:
    GeometryFV3(const eckit::Configuration &, const eckit::mpi::Comm &);
    void print(std::ostream &) const override;

    // FV3-specific accessors
    int npx() const { return npx_; }
    int npy() const { return npy_; }
    int npz() const { return npz_; }
    int tileNum() const { return tileNum_; }

    // Unified access to grid-specific parameters
    eckit::LocalConfiguration gridSpecific() const override;

   private:
    // FV3-specific grid dimensions (cubed-sphere)
    int npx_;     // Number of grid points in x-direction (per tile)
    int npy_;     // Number of grid points in y-direction (per tile)
    int npz_;     // Number of vertical levels (same as numLevels_)
    int tileNum_;  // Tile number for this MPI rank
  };

}  // namespace ijedi

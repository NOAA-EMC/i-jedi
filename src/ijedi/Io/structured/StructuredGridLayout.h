/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "atlas/functionspace.h"
#include "atlas/library/config.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  /// \brief What one MPI task owns of a structured latitude-longitude grid, in the index space
  ///        a netCDF file uses.
  ///
  /// atlas addresses a StructuredColumns function space by (i, j) with j running over
  /// latitudes in the grid's own order (north to south for the Gaussian grids the AI models
  /// use) and i over longitudes. A netCDF file addresses the same points by hyperslab. This
  /// turns the first into the second: the half-open local bounding box [i0, i1) x [j0, j1) is
  /// exactly the block a task has to read, and index(i, j) places a value back into the field.
  ///
  /// For a band-partitioned grid - the usual case - the bounding box is exactly the set of
  /// owned points, so a task reads no more than it needs. For a two-dimensional partitioning
  /// it is a superset, still far smaller than the global field.
  struct StructuredGridLayout
  {
    /// \brief Build from a function space, which must be StructuredColumns over a grid with
    ///        the same number of longitudes on every latitude. Reduced Gaussian grids are
    ///        rejected: they have no rectangular netCDF representation.
    static StructuredGridLayout create(const atlas::FunctionSpace &);

    atlas::idx_t index(atlas::idx_t i, atlas::idx_t j) const { return fs.index(i, j); }

    int boxNx() const { return i1 - i0; }
    int boxNy() const { return j1 - j0; }
    size_t boxSize() const { return static_cast<size_t>(boxNx()) * static_cast<size_t>(boxNy()); }

    atlas::functionspace::StructuredColumns fs;

    int nx = 0;                 ///< global number of longitudes
    int ny = 0;                 ///< global number of latitudes
    std::vector<double> lons;   ///< nx longitudes, atlas order, normalised to [0, 360)
    std::vector<double> lats;   ///< ny latitudes, atlas order

    int i0 = 0, i1 = 0;         ///< local longitude range, half open
    int j0 = 0, j1 = 0;         ///< local latitude range, half open
  };

  // -------------------------------------------------------------------------------------------------

  /// \brief Fold a longitude into [0, 360) so that files written with -180..180 and with
  ///        0..360 compare equal.
  double normaliseLongitude(double lon);

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

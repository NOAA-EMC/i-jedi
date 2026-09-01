/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Io/structured/StructuredGridLayout.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/grid.h"

#include "eckit/exception/Exceptions.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  double normaliseLongitude(double lon)
  {
    double x = std::fmod(lon, 360.0);
    if (x < 0.0) x += 360.0;
    // fmod can land exactly on 360 through rounding; fold it back to 0.
    if (x >= 360.0) x -= 360.0;
    return x;
  }

  // -------------------------------------------------------------------------------------------------

  StructuredGridLayout StructuredGridLayout::create(const atlas::FunctionSpace &functionSpace)
  {
    StructuredGridLayout layout;

    layout.fs = atlas::functionspace::StructuredColumns(functionSpace);
    if (!layout.fs)
    {
      throw eckit::BadValue("StructuredGridLayout: this io backend needs a StructuredColumns "
                            "function space, but the geometry provides '" +
                            functionSpace.type() + "'", Here());
    }

    const atlas::StructuredGrid grid = layout.fs.grid();
    layout.ny = static_cast<int>(grid.ny());
    if (layout.ny <= 0)
    {
      throw eckit::BadValue("StructuredGridLayout: grid has no latitudes", Here());
    }

    layout.nx = static_cast<int>(grid.nx(0));
    for (atlas::idx_t j = 0; j < grid.ny(); ++j)
    {
      if (static_cast<int>(grid.nx(j)) != layout.nx)
      {
        throw eckit::BadValue("StructuredGridLayout: the grid has a different number of "
                              "longitudes on different latitudes (a reduced grid). Such a grid "
                              "has no rectangular netCDF representation; use a regular grid.",
                              Here());
      }
    }

    // Global coordinates, in the grid's own ordering.
    layout.lons.resize(layout.nx);
    for (int i = 0; i < layout.nx; ++i)
    {
      layout.lons[i] = normaliseLongitude(grid.x(i, 0));
    }
    layout.lats.resize(layout.ny);
    for (int j = 0; j < layout.ny; ++j)
    {
      layout.lats[j] = grid.y(j);
    }

    // Local bounding box. j is already contiguous by construction; i is unioned over the
    // owned latitudes, which for a band partitioning spans the full longitude circle.
    layout.j0 = static_cast<int>(layout.fs.j_begin());
    layout.j1 = static_cast<int>(layout.fs.j_end());

    if (layout.j1 <= layout.j0)
    {
      // A task with no owned latitudes still has to take part in collectives, so an empty box
      // is legal; make it explicitly empty rather than inverted.
      layout.j0 = layout.j1 = 0;
      layout.i0 = layout.i1 = 0;
      return layout;
    }

    // Union of the owned longitude ranges. For a band partitioning this is the whole circle;
    // for the two-dimensional partitioning atlas switches to at higher task counts it is a
    // strict subset, and the box is then a superset of what is owned - still far smaller than
    // the global field, and the caller scatters using the per-latitude ranges.
    int iMin = layout.nx;
    int iMax = 0;
    for (int j = layout.j0; j < layout.j1; ++j)
    {
      const int iBegin = static_cast<int>(layout.fs.i_begin(j));
      const int iEnd = static_cast<int>(layout.fs.i_end(j));
      if (iBegin < 0 || iEnd > layout.nx)
      {
        std::stringstream errorMsg;
        errorMsg << "StructuredGridLayout: owned longitude range [" << iBegin << ", " << iEnd
                 << ") on latitude " << j << " falls outside the grid's 0.." << layout.nx
                 << ", which this backend cannot express as a netCDF hyperslab";
        throw eckit::BadValue(errorMsg.str(), Here());
      }
      iMin = std::min(iMin, iBegin);
      iMax = std::max(iMax, iEnd);
    }
    if (iMax <= iMin)
    {
      layout.i0 = layout.i1 = 0;
      return layout;
    }
    layout.i0 = iMin;
    layout.i1 = iMax;

    return layout;
  }

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

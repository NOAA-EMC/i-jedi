/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "ijedi/Geometry/Geometry.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------
/// Maps the owned atlas points of an SQG geometry onto the SQG model's N x N doubly periodic
/// grid (row-major, index iy*N + ix).
///
/// The SQG domain (x, y in [0, L)) is represented in JEDI by a small equatorial lon/lat patch:
/// an N x N regular lon/lat grid with equal spacing in both directions, e.g. an atlas
/// "regional" grid. x increases with longitude and y with latitude; each point's (ix, iy) is
/// recovered from its lon/lat, so the atlas point ordering does not matter. JEDI distances
/// on the patch are the SQG distances scaled by (patch width) / L, so localization radii must
/// be scaled the same way. Periodicity across the patch edges is not represented in JEDI.
struct SQGGrid {
  int N = 0;
  std::vector<int> index;  // per owned point (in geom.ownedNodeIndices() order): iy*N + ix
  std::vector<atlas::idx_t> nodes;  // the owned node indices themselves

  explicit SQGGrid(const Geometry & geom) {
    const atlas::FunctionSpace & fs = geom.functionSpace();
    const auto lonlat = atlas::array::make_view<double, 2>(fs.lonlat());
    nodes = geom.ownedNodeIndices();
    const eckit::mpi::Comm & comm = geom.getComm();

    double lonmin = std::numeric_limits<double>::max(), latmin = lonmin;
    double lonmax = -lonmin, latmax = -lonmin;
    for (const atlas::idx_t n : nodes) {
      lonmin = std::min(lonmin, lonlat(n, 0)); lonmax = std::max(lonmax, lonlat(n, 0));
      latmin = std::min(latmin, lonlat(n, 1)); latmax = std::max(latmax, lonlat(n, 1));
    }
    comm.allReduceInPlace(lonmin, eckit::mpi::min());
    comm.allReduceInPlace(latmin, eckit::mpi::min());
    comm.allReduceInPlace(lonmax, eckit::mpi::max());
    comm.allReduceInPlace(latmax, eckit::mpi::max());

    size_t npts = nodes.size();
    comm.allReduceInPlace(npts, eckit::mpi::sum());
    N = static_cast<int>(std::lround(std::sqrt(static_cast<double>(npts))));
    if (static_cast<size_t>(N) * N != npts || N < 2) {
      throw eckit::BadValue("SQGGrid: geometry has " + std::to_string(npts)
                            + " points, not N*N", Here());
    }
    const double dlon = (lonmax - lonmin) / (N - 1);
    const double dlat = (latmax - latmin) / (N - 1);
    if (std::abs(dlon - dlat) > 1e-6 * dlon) {
      throw eckit::BadValue("SQGGrid: patch spacing must be equal in lon and lat", Here());
    }

    index.resize(nodes.size());
    for (size_t p = 0; p < nodes.size(); ++p) {
      const double fx = (lonlat(nodes[p], 0) - lonmin) / dlon;
      const double fy = (lonlat(nodes[p], 1) - latmin) / dlat;
      const int ix = static_cast<int>(std::lround(fx));
      const int iy = static_cast<int>(std::lround(fy));
      if (std::abs(fx - ix) > 1e-3 || std::abs(fy - iy) > 1e-3) {
        throw eckit::BadValue("SQGGrid: geometry is not a regular lon/lat patch", Here());
      }
      index[p] = iy * N + ix;
    }
  }
};

}  // namespace ijedi

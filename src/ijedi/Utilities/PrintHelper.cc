/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Utilities/PrintHelper.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

#include "atlas/array.h"
#include "atlas/field.h"
#include "oops/util/for_each.h"

namespace ijedi {

std::tuple<double, double, double> fieldMinMaxRMS(const eckit::mpi::Comm & comm,
                                                  const atlas::Field & field,
                                                  const atlas::Field * owned) {
  double localMin          = std::numeric_limits<double>::max();
  double localMax          = std::numeric_limits<double>::lowest();
  double localSumSq        = 0.0;
  atlas::gidx_t localCount = 0;
  const auto fieldView = atlas::array::make_view<double, 2>(field);

  const bool useOwnedMask = owned != nullptr &&
                            owned->rank() == 2 &&
                            owned->shape(0) == field.shape(0) &&
                            owned->shape(1) >= 1;
  // NB: run serially. The reduction below mutates shared accumulators, so the
  // default (parallel) execution pattern would be a data race across OpenMP
  // threads and make the RMS non-reproducible run-to-run (with >1 thread).
  auto accumulateColumn = [&](const atlas::idx_t i) {
    if (useOwnedMask) {
      const auto ownedView = atlas::array::make_view<int, 2>(*owned);
      if (ownedView(i, 0) == 0) {
        return;
      }
    }
    for (atlas::idx_t j = 0; j < field.shape(1); ++j) {
      const double val = fieldView(i, j);
      localMin = std::min(localMin, val);
      localMax = std::max(localMax, val);
      localSumSq += val * val;
      ++localCount;
    }
  };

  if (!util::details::hasContiguousOwnedPoints(field)) {
    const auto ghost = atlas::array::make_view<int, 1>(field.functionspace().ghost());
    for (atlas::idx_t i = 0; i < field.shape(0); ++i) {
      if (ghost(i) == 0) {
        accumulateColumn(i);
      }
    }
  } else {
    const auto indexSpace = util::make_index_space_2d(util::IndexRange::exclude_halo, field);
    for (atlas::idx_t i = indexSpace.i.start; i < indexSpace.i.end; ++i) {
      accumulateColumn(i);
    }
  }
  double globalMin          = localMin;
  double globalMax          = localMax;
  double globalSumSq        = localSumSq;
  atlas::gidx_t globalCount = localCount;
  comm.allReduceInPlace(globalMin, eckit::mpi::min());
  comm.allReduceInPlace(globalMax, eckit::mpi::max());
  comm.allReduceInPlace(globalSumSq, eckit::mpi::sum());
  comm.allReduceInPlace(globalCount, eckit::mpi::sum());
  if (globalCount == 0) {
    return std::make_tuple(0.0, 0.0, 0.0);
  }
  double rms = (globalCount > 0) ? std::sqrt(globalSumSq / globalCount) : 0.0;
  return std::make_tuple(globalMin, globalMax, rms);
}

}  // namespace ijedi

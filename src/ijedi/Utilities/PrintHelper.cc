/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Utilities/PrintHelper.h"

#include <cmath>
#include <limits>
#include <tuple>

#include "atlas/field.h"
#include "oops/util/for_each.h"

namespace ijedi {

std::tuple<double, double, double> fieldMinMaxRMS(const eckit::mpi::Comm & comm,
                                                  const atlas::Field & field) {
  double localMin          = std::numeric_limits<double>::max();
  double localMax          = std::numeric_limits<double>::lowest();
  double localSumSq        = 0.0;
  atlas::gidx_t localCount = 0;
  util::for_each_value(
      util::IndexRange::exclude_halo,
      [&](const double val) {
        localMin = std::min(localMin, val);
        localMax = std::max(localMax, val);
        localSumSq += val * val;
        ++localCount;
      },
      field);
  double globalMin          = localMin;
  double globalMax          = localMax;
  double globalSumSq        = localSumSq;
  atlas::gidx_t globalCount = localCount;
  comm.allReduceInPlace(globalMin, eckit::mpi::min());
  comm.allReduceInPlace(globalMax, eckit::mpi::max());
  comm.allReduceInPlace(globalSumSq, eckit::mpi::sum());
  comm.allReduceInPlace(globalCount, eckit::mpi::sum());
  double rms = (globalCount > 0) ? std::sqrt(globalSumSq / globalCount) : 0.0;
  return std::make_tuple(globalMin, globalMax, rms);
}

}  // namespace ijedi
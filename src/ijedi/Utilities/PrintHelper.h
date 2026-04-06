/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <tuple>

namespace atlas {
class Field;
}  // namespace atlas

namespace eckit::mpi {
class Comm;
}  // namespace eckit::mpi

namespace ijedi {

std::tuple<double, double, double> fieldMinMaxRMS(const eckit::mpi::Comm &, const atlas::Field &);

}  // namespace ijedi
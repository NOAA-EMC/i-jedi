/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Increment/Increment.h"

#include <algorithm>
#include <cmath>
#include <iomanip>

#include "atlas/field.h"
#include "eckit/config/Configuration.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Utilities/PrintHelper.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/for_each.h"

namespace ijedi {

Increment::Increment(const Geometry & geom, const oops::Variables & vars,
                     const util::DateTime & time)
    : mist::base::Increment(geom, vars, time), geom_(geom) {}

Increment::Increment(const Geometry & geom, const Increment & other, const bool ad)
    : mist::base::Increment(geom, other, ad), geom_(geom) {}

Increment::Increment(const Increment & other, const bool copy)
    : mist::base::Increment(other, copy), geom_(other.geom_) {}

Increment::~Increment() = default;

Increment & Increment::operator=(const Increment & rhs) {
  mist::base::Increment::operator=(rhs);
  return *this;
}

void Increment::read(const eckit::Configuration & config) {
  oops::Log::trace() << "ijedi::Increment::read starting" << std::endl;

  // TODO(someone) implement read

  oops::Log::trace() << "ijedi::Increment::read done" << std::endl;
}

void Increment::write(const eckit::Configuration & config) const {
  oops::Log::trace() << "ijedi::Increment::write starting" << std::endl;

  // TODO(someone) implement write

  oops::Log::trace() << "ijedi::Increment::write done" << std::endl;
}

void Increment::print(std::ostream & os) const {
  os << std::endl
     << "  Valid time: " << this->validTime() << std::endl
     << ", nFields = " << this->variables().size();

  const auto & comm = geom_.comm();
  const auto & fs   = this->fieldSet();
  int fieldIndex    = 0;
  for (const auto & var : this->variables()) {
    const atlas::Field & field            = fs.field(var.name());
    const auto[globalMin, globalMax, rms] = fieldMinMaxRMS(comm, field);
    ++fieldIndex;
    os << std::endl
       << "Fld=" << fieldIndex << std::scientific << std::setprecision(16) << "  Min=" << globalMin
       << ", Max=" << globalMax << ", RMS=" << rms << " : " << var.name();
  }
}

}  // namespace ijedi

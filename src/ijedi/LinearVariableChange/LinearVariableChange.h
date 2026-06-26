/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "mist/base/LinearVariableChange.h"

namespace ijedi {

class Geometry;

class LinearVariableChange : public mist::LinearVariableChange {
 public:
  LinearVariableChange(const Geometry &, const eckit::Configuration &);
};

}  // namespace ijedi

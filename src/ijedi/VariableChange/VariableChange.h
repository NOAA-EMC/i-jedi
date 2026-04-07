/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "mist/base/VariableChange.h"

namespace ijedi {

class Geometry;

class VariableChange : public mist::base::VariableChange {
 public:
  VariableChange(const eckit::Configuration &, const Geometry & geometry);
};

}  // namespace ijedi

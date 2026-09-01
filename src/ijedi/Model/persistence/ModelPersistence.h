/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "ijedi/Model/base/ModelBase.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace mist {
class ModelAuxControl;
}  // namespace mist

namespace ijedi {

class Geometry;
class State;

// -------------------------------------------------------------------------------------------------

/// \brief Persistence: the state is carried forward unchanged, only its valid time advances.
///
/// Useful as a null model for testing the forecast machinery, and as the trivial baseline a
/// real forecast has to beat.
class ModelPersistence : public ModelBase {
 public:
  static const std::string classname() { return "ijedi::ModelPersistence"; }

  ModelPersistence(const Geometry & geom, const eckit::Configuration & config)
    : ModelBase(geom, config) {}
  ~ModelPersistence() = default;

  void initialize(State &) const override {}
  void step(State &, const mist::ModelAuxControl &) const override;
  void finalize(State &) const override {}

 private:
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "ijedi/Model/python/PyModelBridge.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace ijedi {

// -------------------------------------------------------------------------------------------------

/// \brief A bridge backend with no Python behind it: the declaration comes from the yaml and
///        the model holds the fields it was given.
///
/// This exists so that everything on the C++ side of the boundary - the declaration contract,
/// the geometry and variable checks, the gather and scatter, the forecast lifecycle - can be
/// run and tested before any framework backend exists, and so that a backend author has a
/// worked example of the interface. Its dynamics are persistence, which also makes it a
/// useful null model when something else is under test.
class StubPyModelBridge : public PyModelBridge {
 public:
  static const std::string classname() { return "ijedi::StubPyModelBridge"; }

  explicit StubPyModelBridge(const eckit::Configuration & config);
  ~StubPyModelBridge() = default;

  const Declaration & declaration() const override { return declaration_; }

  void encode(const Fields & inputs, const util::DateTime & validTime) override;
  void advance(int steps) override;
  void decode(Fields & outputs) const override;
  void reset() override;

 private:
  Declaration declaration_;
  Fields state_;
  bool encoded_ = false;
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

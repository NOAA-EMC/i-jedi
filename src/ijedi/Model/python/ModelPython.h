/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <ostream>
#include <string>

#include "ijedi/Model/base/ModelBase.h"
#include "ijedi/Model/python/PyModelBridge.h"

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

/// \brief Any forecast model written in Python, driven through PyModelBridge.
///
/// There is deliberately no C++ class per model. NeuralGCM, AIGFS and a PyTorch physical model
/// differ in ways that are entirely model knowledge - which variables, on which levels, under
/// which names, with which internal timestep - and that knowledge belongs in the Python
/// adapter, which reports it through PyModelBridge::Declaration. This class checks the
/// declaration against the geometry and does the marshalling. Adding a model is a Python-only
/// change: write an adapter, name it in the yaml.
///
/// Models run on their own native grid; i-jedi interpolates nothing here. Configure the
/// geometry to match what the adapter declares and the checks below will confirm it.
///
/// The lifecycle maps onto the adapter's:
///   constructor   ->  load the model, read its declaration, check it against the geometry
///   initialize()  ->  encode    build the internal state from the current State
///   step()        ->  advance   run the internal steps that make up one oops step, then decode
///   finalize()    ->  reset     drop the internal state, keep the model loaded
/// so "time step" in the yaml is the interval at which oops sees the state, not the model's
/// own internal timestep, which the adapter reports and which is usually far shorter.
class ModelPython : public ModelBase {
 public:
  static const std::string classname() { return "ijedi::ModelPython"; }

  ModelPython(const Geometry &, const eckit::Configuration &);
  ~ModelPython() = default;

  void initialize(State &) const override;
  void step(State &, const mist::ModelAuxControl &) const override;
  void finalize(State &) const override;

  /// \brief JEDI long name -> adapter variable name, as resolved from the adapter's
  ///        declaration and the yaml's "field name map" override. Exposed for testing.
  const std::map<std::string, std::string> & fieldNameMap() const { return fieldNameMap_; }

 private:
  void print(std::ostream &) const override;

  /// Check the adapter's declared levels and variables against the geometry, so that a
  /// mismatch is a startup error rather than a strange forecast.
  void checkDeclaration(const PyModelBridge::Declaration &);

  const Geometry & geom_;

  std::string adapter_;
  std::string backend_;
  std::map<std::string, std::string> fieldNameMap_;
  std::map<std::string, double> fieldScaling_;

  // The model outlives a forecast leg; the state it encodes does not.
  std::unique_ptr<PyModelBridge> bridge_;
  int innerSteps_ = 0;
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/Duration.h"

#include "ijedi/Model/ModelBase.h"
#include "ijedi/Model/sqg/SQGGrid.h"

class SQG;  // Model/sqg/sqg.hpp

namespace ijedi {

class Geometry;
class State;

// -------------------------------------------------------------------------------------------------
/// Surface quasi-geostrophic (SQG) turbulence model (Model/sqg/sqg.hpp) as an ijedi forecast
/// model, `model: name: SQG`.
///
/// State: "air_potential_temperature" (K) on 2 levels, level 0 the lower boundary (z = 0) and
/// level 1 the lid (z = H) -- the SQG boundary PV scaled by f*theta0/g. The geometry is an
/// N x N equatorial lon/lat patch (see SQGGrid). Each step converts the state to SQG PV,
/// runs (time step)/dt SQG steps and converts back; the conversion round-trips SQG's float
/// PV exactly. Runs on a single MPI task (the model is spectral over the whole domain).
///
/// YAML (SQG parameters, defaults as in sqgturb-cpp/main.cpp):
///   name: SQG
///   time step: PT6H           # JEDI step; must be a multiple of dt
///   dt: 1800                  # SQG time step (s)
///   L: 20.e6, H: 10.e3, U: 20, f: 1.e-4, nsq: 1.e-4, r: 0, tdiab: 864000,
///   diff_order: 8, diff_efold: 86400, theta0: 300, g: 9.8
///   nthreads: 1               # OpenMP/FFT threads (0: OMP_NUM_THREADS)
class ModelSQG : public ModelBase {
 public:
  static const std::string classname() { return "ijedi::ModelSQG"; }

  ModelSQG(const Geometry &, const eckit::Configuration &, const util::Duration &);
  ~ModelSQG() override;

  void step(State &) const override;
  const oops::Variables & variables() const override { return vars_; }
  void print(std::ostream &) const override;

  static constexpr const char * varname = "air_potential_temperature";

 private:
  SQGGrid grid_;
  int nsteps_;           // SQG steps per JEDI step
  float scale_;          // f*theta0/g: SQG PV (m/s) -> potential temperature (K)
  float dt_;
  int nthreads_;
  oops::Variables vars_;
  mutable std::unique_ptr<SQG> sqg_;
};

}  // namespace ijedi

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
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

#include "ijedi/LinearModel/LinearModelBase.h"
#include "ijedi/Model/sqg/SQGGrid.h"

class SQGTL;  // Model/sqg/sqg_tl.hpp

namespace ijedi {

class Geometry;
class Increment;
class State;

// -------------------------------------------------------------------------------------------------
/// Tangent linear and adjoint of the SQG model (ModelSQG), `linear model: name: SQG`, e.g. for
/// 4D-Var. Takes the same SQG parameters as ModelSQG (they must match the nonlinear model)
/// plus `time step`. Each TL step linearizes ModelSQG's step about the trajectory state saved
/// at its start time; TL and adjoint run in double precision (Model/sqg/sqg_tl.hpp) and the
/// adjoint is exact for the Euclidean inner product of increments. Single MPI task.
class LinearModelSQG : public LinearModelBase {
 public:
  static const std::string classname() { return "ijedi::LinearModelSQG"; }

  LinearModelSQG(const Geometry &, const eckit::Configuration &, const util::Duration &);
  ~LinearModelSQG() override;

  void setTrajectory(const State &) override;
  void stepTL(Increment &) const override;
  void stepAD(Increment &) const override;
  void print(std::ostream &) const override;

 private:
  const std::vector<double> & trajectory(const util::DateTime &) const;
  std::vector<double> gather(const Increment &) const;
  void scatter(const std::vector<double> &, Increment &) const;

  SQGGrid grid_;
  util::Duration tstep_;
  int nsteps_;
  float scale_;
  std::unique_ptr<SQGTL> tl_;
  std::map<util::DateTime, std::vector<double>> traj_;  // grid PV (m/s) at each step start
};

}  // namespace ijedi

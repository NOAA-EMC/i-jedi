/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/sqg/LinearModelSQG.h"

#include <cmath>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/Model/sqg/ModelSQG.h"
#include "ijedi/Model/sqg/sqg_tl.hpp"
#include "ijedi/State/State.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------
static LinearModelMaker<LinearModelSQG> makerLinearModelSQG_("SQG");

// -------------------------------------------------------------------------------------------------
LinearModelSQG::LinearModelSQG(const Geometry & geom, const eckit::Configuration & conf,
                               const util::Duration & tstep)
  : grid_(geom), tstep_(tstep) {
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  if (geom.getComm().size() != 1) {
    throw eckit::NotImplemented("LinearModelSQG runs on a single MPI task", Here());
  }
  if (geom.numLevels() != 2) {
    throw eckit::BadValue("LinearModelSQG: geometry must have 2 levels (z=0 and z=H)", Here());
  }
  // Parameters and defaults as ModelSQG (read as float there, so convert the same way)
  const float dt = conf.getFloat("dt", 1800.f);
  const double nst = static_cast<double>(tstep.toSeconds()) / dt;
  nsteps_ = static_cast<int>(std::lround(nst));
  if (nsteps_ < 1 || std::abs(nst - nsteps_) > 1e-9) {
    throw eckit::BadValue("LinearModelSQG: time step must be a multiple of dt", Here());
  }
  const float f = conf.getFloat("f", 1.e-4f), theta0 = conf.getFloat("theta0", 300.f);
  const float g = conf.getFloat("g", 9.8f);
  scale_ = f * theta0 / g;
  tl_ = std::make_unique<SQGTL>(
      grid_.N, f, conf.getFloat("nsq", 1.e-4f), conf.getFloat("L", 20.e6f),
      conf.getFloat("H", 10.e3f), conf.getFloat("U", 20.f), conf.getFloat("r", 0.f),
      conf.getFloat("tdiab", 10.f * 86400.f), conf.getInt("diff_order", 8),
      conf.getFloat("diff_efold", 86400.f), dt);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
LinearModelSQG::~LinearModelSQG() = default;

// -------------------------------------------------------------------------------------------------
void LinearModelSQG::setTrajectory(const State & xx) {
  // The SQG PV the nonlinear step starts from (ModelSQG::step's float conversion).
  const int NN = grid_.N * grid_.N;
  const auto view = atlas::array::make_view<double, 2>(xx.fieldSet().field(ModelSQG::varname));
  const double scale = scale_;
  std::vector<double> pv(2 * NN);
  for (size_t p = 0; p < grid_.nodes.size(); ++p) {
    for (int k = 0; k < 2; ++k) {
      pv[k * NN + grid_.index[p]] = static_cast<float>(view(grid_.nodes[p], k) / scale);
    }
  }
  traj_[xx.validTime()] = std::move(pv);
}

// -------------------------------------------------------------------------------------------------
const std::vector<double> & LinearModelSQG::trajectory(const util::DateTime & time) const {
  const auto it = traj_.find(time);
  if (it == traj_.end()) {
    throw eckit::BadValue("LinearModelSQG: no trajectory at " + time.toString(), Here());
  }
  return it->second;
}

// -------------------------------------------------------------------------------------------------
// Increments stay in K: the PV <-> K scale factor cancels in a linear map.
std::vector<double> LinearModelSQG::gather(const Increment & dx) const {
  const int NN = grid_.N * grid_.N;
  const auto view = atlas::array::make_view<double, 2>(dx.fieldSet().field(ModelSQG::varname));
  std::vector<double> v(2 * NN);
  for (size_t p = 0; p < grid_.nodes.size(); ++p) {
    for (int k = 0; k < 2; ++k) v[k * NN + grid_.index[p]] = view(grid_.nodes[p], k);
  }
  return v;
}

void LinearModelSQG::scatter(const std::vector<double> & v, Increment & dx) const {
  const int NN = grid_.N * grid_.N;
  auto view = atlas::array::make_view<double, 2>(dx.fieldSet().field(ModelSQG::varname));
  for (size_t p = 0; p < grid_.nodes.size(); ++p) {
    for (int k = 0; k < 2; ++k) view(grid_.nodes[p], k) = v[k * NN + grid_.index[p]];
  }
}

// -------------------------------------------------------------------------------------------------
void LinearModelSQG::stepTL(Increment & dx) const {
  oops::Log::trace() << classname() << " stepTL starting" << std::endl;
  scatter(tl_->tl(trajectory(dx.validTime()), gather(dx), nsteps_), dx);
  dx.updateTime(tstep_);
  oops::Log::trace() << classname() << " stepTL done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
void LinearModelSQG::stepAD(Increment & dx) const {
  oops::Log::trace() << classname() << " stepAD starting" << std::endl;
  const util::DateTime start = dx.validTime() - tstep_;
  scatter(tl_->ad(trajectory(start), gather(dx), nsteps_), dx);
  dx.updateTime(-tstep_);
  oops::Log::trace() << classname() << " stepAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
void LinearModelSQG::print(std::ostream & os) const {
  os << "SQG tangent linear and adjoint model: N=" << grid_.N << ", " << nsteps_
     << " SQG steps per model step of " << tstep_ << ", trajectory at " << traj_.size()
     << " time(s)";
}

}  // namespace ijedi

/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/sqg/ModelSQG.h"

#include <cmath>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Model/sqg/sqg.hpp"
#include "ijedi/State/State.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------
static ModelMaker<ModelSQG> makerModelSQG_("SQG");

// -------------------------------------------------------------------------------------------------
ModelSQG::ModelSQG(const Geometry & geom, const eckit::Configuration & conf,
                   const util::Duration & tstep)
  : grid_(geom), vars_(std::vector<std::string>{varname}) {
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  if (geom.getComm().size() != 1) {
    throw eckit::NotImplemented("ModelSQG runs on a single MPI task", Here());
  }
  if (geom.numLevels() != 2) {
    throw eckit::BadValue("ModelSQG: geometry must have 2 levels (z=0 and z=H)", Here());
  }

  dt_ = conf.getFloat("dt", 1800.f);
  const double nst = static_cast<double>(tstep.toSeconds()) / dt_;
  nsteps_ = static_cast<int>(std::lround(nst));
  if (nsteps_ < 1 || std::abs(nst - nsteps_) > 1e-9) {
    throw eckit::BadValue("ModelSQG: time step must be a multiple of dt", Here());
  }
  const float f = conf.getFloat("f", 1.e-4f), theta0 = conf.getFloat("theta0", 300.f);
  const float g = conf.getFloat("g", 9.8f);
  scale_ = f * theta0 / g;  // in float, as sqgturb-cpp/main.cpp and IoSQG
  nthreads_ = conf.getInt("nthreads", 1);

  const int N = grid_.N;
  const std::vector<float> pv0(2 * N * N, 0.f);  // replaced by the state at every step
  sqg_ = std::make_unique<SQG>(
      pv0, N, f, conf.getFloat("nsq", 1.e-4f), conf.getFloat("L", 20.e6f),
      conf.getFloat("H", 10.e3f), conf.getFloat("U", 20.f), conf.getFloat("r", 0.f),
      conf.getFloat("tdiab", 10.f * 86400.f), conf.getInt("diff_order", 8),
      conf.getFloat("diff_efold", 86400.f), theta0, g, dt_, 0.0, nthreads_);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
ModelSQG::~ModelSQG() = default;

// -------------------------------------------------------------------------------------------------
void ModelSQG::step(State & xx) const {
  oops::Log::trace() << classname() << " step starting" << std::endl;
  const int NN = grid_.N * grid_.N;
  auto view = atlas::array::make_view<double, 2>(xx.fieldSet().field(varname));
  const double scale = scale_;

  std::vector<float> pv(2 * NN);
  for (size_t p = 0; p < grid_.nodes.size(); ++p) {
    for (int k = 0; k < 2; ++k) {
      pv[k * NN + grid_.index[p]] = static_cast<float>(view(grid_.nodes[p], k) / scale);
    }
  }
  const std::vector<float> out = sqg_->advance(nsteps_, &pv);
  for (size_t p = 0; p < grid_.nodes.size(); ++p) {
    for (int k = 0; k < 2; ++k) {
      view(grid_.nodes[p], k) = static_cast<double>(out[k * NN + grid_.index[p]]) * scale;
    }
  }
  oops::Log::trace() << classname() << " step done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
void ModelSQG::print(std::ostream & os) const {
  os << "SQG model: N=" << grid_.N << ", dt=" << dt_ << " s, " << nsteps_
     << " SQG steps per model step, " << nthreads_ << " thread(s), scale to K " << scale_;
}

}  // namespace ijedi

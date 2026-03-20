// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/State.h"
#include "ijedi/State/base/StateBase.h"

namespace ijedi
{

  // ---------------------------------------------------------------------------
  State::State(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & time)
  {
    oops::Log::trace() << "State::State starting" << std::endl;
    stateImpl_ = StateBase::create(geom, vars, time);
    oops::Log::trace() << "State::State done" << std::endl;
  }

  // ---------------------------------------------------------------------------
  State::State(const Geometry & geom, const eckit::Configuration & conf)
  {
    oops::Log::trace() << "State::State (config) starting" << std::endl;
    stateImpl_ = StateBase::create(geom, conf);
    oops::Log::trace() << "State::State (config) done" << std::endl;
  }

  // ---------------------------------------------------------------------------
  State::State(const Geometry & geom, const State & other)
  {
    oops::Log::trace() << "State::State (change geom) starting" << std::endl;
    stateImpl_ = StateBase::create(geom, *other.stateImpl_);
    oops::Log::trace() << "State::State (change geom) done" << std::endl;
  }

  // ---------------------------------------------------------------------------
  State::State(const oops::Variables & vars, const State & other)
  {
    oops::Log::trace() << "State::State (change vars) starting" << std::endl;
    stateImpl_ = StateBase::create(vars, *other.stateImpl_);
    oops::Log::trace() << "State::State (change vars) done" << std::endl;
  }

  // ---------------------------------------------------------------------------
  State::State(const State & other)
  {
    oops::Log::trace() << "State::State (copy) starting" << std::endl;
    stateImpl_ = other.stateImpl_;
    oops::Log::trace() << "State::State (copy) done" << std::endl;
  }

  // ---------------------------------------------------------------------------
  State::~State() {}

  // ---------------------------------------------------------------------------
  State & State::operator=(const State & rhs)
  {
    stateImpl_->copyDataFrom(*rhs.stateImpl_);
    return *this;
  }

  // ---------------------------------------------------------------------------
  const util::DateTime & State::validTime() const
  { return stateImpl_->validTime(); }

  void State::updateTime(const util::Duration & dt)
  { stateImpl_->updateTime(dt); }

  // ---------------------------------------------------------------------------
  void State::zero()                              { stateImpl_->zero(); }
  void State::accumul(const double & w, const State & x)
  { stateImpl_->accumul(w, *x.stateImpl_); }
  double State::norm() const                      { return stateImpl_->norm(); }

  // ---------------------------------------------------------------------------
  void State::read(const eckit::Configuration & conf)
  { stateImpl_->read(conf); }
  void State::write(const eckit::Configuration & conf) const
  { stateImpl_->write(conf); }

  // ---------------------------------------------------------------------------
  size_t State::serialSize() const
  { return stateImpl_->serialSize(); }

  void State::serialize(std::vector<double> & vect) const
  { stateImpl_->serialize(vect); }

  void State::deserialize(const std::vector<double> & vect, size_t & current)
  { stateImpl_->deserialize(vect, current); }

  // ---------------------------------------------------------------------------
  void State::toFieldSet(atlas::FieldSet & fset) const
  { stateImpl_->toFieldSet(fset); }

  void State::fromFieldSet(const atlas::FieldSet & fset)
  { stateImpl_->fromFieldSet(fset); }

  // ---------------------------------------------------------------------------
  const oops::Variables & State::variables() const
  { return stateImpl_->variables(); }

  // ---------------------------------------------------------------------------
  void State::print(std::ostream & os) const
  { stateImpl_->print(os); }

}  // namespace ijedi

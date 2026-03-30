// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/Increment/base/IncrementBase.h"
#include "ijedi/State/State.h"

namespace ijedi
{

  // ---------------------------------------------------------------------------
  Increment::Increment(const Geometry & geom, const oops::Variables & vars,
                       const util::DateTime & time)
  {
    oops::Log::trace() << "Increment::Increment starting" << std::endl;
    impl_ = IncrementBase::create(geom, vars, time);
    oops::Log::trace() << "Increment::Increment done" << std::endl;
  }

  // ---------------------------------------------------------------------------
  Increment::Increment(const Geometry & geom, const Increment & other)
  {
    oops::Log::trace() << "Increment::Increment (change geom) starting"
                       << std::endl;
    impl_ = IncrementBase::create(geom, *other.impl_);
    oops::Log::trace() << "Increment::Increment (change geom) done"
                       << std::endl;
  }

  // ---------------------------------------------------------------------------
  Increment::Increment(const Increment & other, bool copy)
  {
    oops::Log::trace() << "Increment::Increment (copy) starting" << std::endl;
    impl_ = IncrementBase::create(*other.impl_, copy);
    oops::Log::trace() << "Increment::Increment (copy) done" << std::endl;
  }

  // ---------------------------------------------------------------------------
  Increment::~Increment() {}

  // ---------------------------------------------------------------------------
  Increment & Increment::operator=(const Increment & rhs)
  {
    impl_->copyDataFrom(*rhs.impl_);
    return *this;
  }

  // ---------------------------------------------------------------------------
  const util::DateTime & Increment::validTime() const
  { return impl_->validTime(); }

  void Increment::updateTime(const util::Duration & dt)
  { impl_->updateTime(dt); }

  // ---------------------------------------------------------------------------
  void Increment::diff(const State & x1, const State & x2)
  { impl_->diff(x1.impl(), x2.impl()); }

  void Increment::zero()                    { impl_->zero(); }
  void Increment::zero(const util::DateTime & dt) { impl_->zero(dt); }
  void Increment::ones()                    { impl_->ones(); }
  void Increment::sqrt()                    { impl_->sqrt(); }
  void Increment::random()                  { impl_->random(); }

  void Increment::dirac(const eckit::Configuration & conf)
  { impl_->dirac(conf); }

  void Increment::accumul(const double & w, const State & x)
  { impl_->accumul(w, x.impl()); }

  // ---------------------------------------------------------------------------
  Increment & Increment::operator+=(const Increment & rhs)
  { *impl_ += *rhs.impl_; return *this; }

  Increment & Increment::operator-=(const Increment & rhs)
  { *impl_ -= *rhs.impl_; return *this; }

  Increment & Increment::operator*=(const double & scalar)
  { *impl_ *= scalar; return *this; }

  void Increment::axpy(const double & w, const Increment & rhs, const bool /*check*/)
  { impl_->axpy(w, *rhs.impl_); }

  double Increment::dot_product_with(const Increment & other) const
  { return impl_->dot_product_with(*other.impl_); }

  void Increment::schur_product_with(const Increment & other)
  { impl_->schur_product_with(*other.impl_); }

  // ---------------------------------------------------------------------------
  void Increment::read(const eckit::Configuration & conf)
  { impl_->read(conf); }

  void Increment::write(const eckit::Configuration & conf) const
  { impl_->write(conf); }

  double Increment::norm() const
  { return impl_->norm(); }

  // ---------------------------------------------------------------------------
  size_t Increment::serialSize() const
  { return impl_->serialSize(); }

  void Increment::serialize(std::vector<double> & vect) const
  { impl_->serialize(vect); }

  void Increment::deserialize(const std::vector<double> & vect, size_t & current)
  { impl_->deserialize(vect, current); }

  // ---------------------------------------------------------------------------
  void Increment::toFieldSet(atlas::FieldSet & fset) const
  { impl_->toFieldSet(fset); }

  void Increment::fromFieldSet(const atlas::FieldSet & fset)
  { impl_->fromFieldSet(fset); }

  // ---------------------------------------------------------------------------
  const oops::Variables & Increment::variables() const
  { return impl_->variables(); }

  // ---------------------------------------------------------------------------
  void Increment::print(std::ostream & os) const
  { impl_->print(os); }

}  // namespace ijedi

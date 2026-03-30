// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <string>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Increment/base/IncrementBase.h"
#include "ijedi/Increment/mom6/IncrementMOM6.h"
#include "ijedi/State/base/StateBase.h"

namespace ijedi
{

  // ---------------------------------------------------------------------------
  // Factory helpers
  // ---------------------------------------------------------------------------
  namespace
  {
    std::string geomType(const Geometry & geom)
    {
      return geom.gridSpecific().getString("grid_type");
    }
  }

  std::shared_ptr<IncrementBase> IncrementBase::create(
      const Geometry & geom,
      const oops::Variables & vars,
      const util::DateTime & time)
  {
    const std::string type = geomType(geom);
    if (type == "mom6")
      return std::make_shared<IncrementMOM6>(geom, vars, time);
    throw eckit::BadValue("Unsupported geometry type for Increment: " + type,
                          Here());
  }

  std::shared_ptr<IncrementBase> IncrementBase::create(
      const Geometry & geom,
      const IncrementBase & other)
  {
    const std::string type = geomType(geom);
    if (type == "mom6")
      return std::make_shared<IncrementMOM6>(geom, other);
    throw eckit::BadValue("Unsupported geometry type for Increment: " + type,
                          Here());
  }

  std::shared_ptr<IncrementBase> IncrementBase::create(
      const IncrementBase & other, bool copy)
  {
    if (const auto * dx = dynamic_cast<const IncrementMOM6 *>(&other))
      return std::make_shared<IncrementMOM6>(*dx, copy);
    throw eckit::BadValue(
        "Unsupported IncrementBase subtype in copy constructor", Here());
  }

  // ---------------------------------------------------------------------------
  // Private helper — allocate one zero-filled field per variable
  // ---------------------------------------------------------------------------
  void IncrementBase::allocateFields()
  {
    for (int iv = 0; iv < static_cast<int>(vars_.size()); ++iv)
    {
      const std::string & name = vars_[iv].name();
      atlas::Field f = functionSpace_.createField<double>(
          atlas::option::name(name) | atlas::option::levels(numLevels_));
      auto view = atlas::array::make_view<double, 2>(f);
      view.assign(0.0);
      fields_.add(f);
    }
  }

  // ---------------------------------------------------------------------------
  // Protected init helpers
  // ---------------------------------------------------------------------------
  void IncrementBase::initFromGeomVarsTime(const Geometry & geom,
                                           const oops::Variables & vars,
                                           const util::DateTime & time)
  {
    functionSpace_ = geom.functionSpace();
    numLevels_     = geom.numLevels();
    vars_          = vars;
    time_          = time;
    allocateFields();
  }

  void IncrementBase::initFromGeomIncrement(const Geometry & geom,
                                            const IncrementBase & other)
  {
    functionSpace_ = geom.functionSpace();
    numLevels_     = geom.numLevels();
    vars_          = other.vars_;
    time_          = other.time_;
    allocateFields();

    for (int iv = 0; iv < static_cast<int>(vars_.size()); ++iv)
    {
      const std::string & name = vars_[iv].name();
      if (!other.fields_.has(name)) continue;
      auto       dst = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src = atlas::array::make_view<double, 2>(
                           other.fields_[name]);
      const int npts = std::min(fields_[iv].shape(0),
                                other.fields_[name].shape(0));
      const int nlev = std::min(fields_[iv].shape(1),
                                other.fields_[name].shape(1));
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) = src(n, l);
    }
  }

  void IncrementBase::initCopy(const IncrementBase & other, bool copy)
  {
    functionSpace_ = other.functionSpace_;
    numLevels_     = other.numLevels_;
    vars_          = other.vars_;
    time_          = other.time_;
    allocateFields();
    if (copy) copyDataFrom(other);
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::copyDataFrom(const IncrementBase & other)
  {
    vars_ = other.vars_;
    time_ = other.time_;
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      auto       dst = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src = atlas::array::make_view<double, 2>(
                           other.fields_[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) = src(n, l);
    }
  }

  // ---------------------------------------------------------------------------
  // Generic Atlas-based operations
  // ---------------------------------------------------------------------------
  void IncrementBase::zero()
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      auto view = atlas::array::make_view<double, 2>(fields_[iv]);
      view.assign(0.0);
    }
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::ones()
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      auto view = atlas::array::make_view<double, 2>(fields_[iv]);
      view.assign(1.0);
    }
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::sqrt()
  {
    const auto ghost = atlas::array::make_view<int32_t, 1>(
                           functionSpace_.ghost());
    // First pass: check for negative values (throws before modifying any data)
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const auto view = atlas::array::make_view<double, 2>(fields_[iv]);
      const int npts  = fields_[iv].shape(0);
      const int nlev  = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
      {
        if (ghost(n)) continue;
        for (int l = 0; l < nlev; ++l)
          if (view(n, l) < 0.0)
            throw eckit::BadValue(
                "IncrementBase::sqrt — negative value in increment field "
                + vars_[iv].name(), Here());
      }
    }
    // Second pass: apply sqrt
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      auto view      = atlas::array::make_view<double, 2>(fields_[iv]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
      {
        if (ghost(n)) continue;
        for (int l = 0; l < nlev; ++l)
          view(n, l) = std::sqrt(view(n, l));
      }
    }
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::random()
  {
    const auto ghost = atlas::array::make_view<int32_t, 1>(
                           functionSpace_.ghost());
    std::mt19937 gen(std::random_device{}());
    std::normal_distribution<double> dist(0.0, 1.0);

    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      auto view      = atlas::array::make_view<double, 2>(fields_[iv]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
      {
        if (ghost(n)) { continue; }
        for (int l = 0; l < nlev; ++l)
          view(n, l) = dist(gen);
      }
    }
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::dirac(const eckit::Configuration & conf)
  {
    zero();

    const int targetIdx = conf.getInt("indices");
    const int targetLev = conf.getInt("level", 0);
    const std::string var = conf.getString("variable", vars_[0].name());

    if (!fields_.has(var))
      throw eckit::BadValue("dirac: unknown variable " + var, Here());

    auto view = atlas::array::make_view<double, 2>(fields_[var]);
    const int npts = fields_[var].shape(0);
    const int nlev = fields_[var].shape(1);

    ASSERT_MSG(targetIdx < npts,
        "dirac: indices " + std::to_string(targetIdx) +
        " out of range [0, " + std::to_string(npts - 1) + "]");
    ASSERT_MSG(targetLev < nlev,
        "dirac: level " + std::to_string(targetLev) +
        " out of range [0, " + std::to_string(nlev - 1) + "]");

    view(targetIdx, targetLev) = 1.0;
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::diff(const StateBase & x1, const StateBase & x2)
  {
    time_ = x1.validTime();
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      auto       dx  = atlas::array::make_view<double, 2>(fields_[name]);
      const auto v1  = atlas::array::make_view<double, 2>(x1.fields()[name]);
      const auto v2  = atlas::array::make_view<double, 2>(x2.fields()[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dx(n, l) = v1(n, l) - v2(n, l);
    }
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::accumul(const double & w, const StateBase & x)
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      auto       dst = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src = atlas::array::make_view<double, 2>(x.fields()[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) += w * src(n, l);
    }
  }

  // ---------------------------------------------------------------------------
  IncrementBase & IncrementBase::operator+=(const IncrementBase & rhs)
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      auto       dst = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src = atlas::array::make_view<double, 2>(rhs.fields_[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) += src(n, l);
    }
    return *this;
  }

  // ---------------------------------------------------------------------------
  IncrementBase & IncrementBase::operator-=(const IncrementBase & rhs)
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      auto       dst = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src = atlas::array::make_view<double, 2>(rhs.fields_[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) -= src(n, l);
    }
    return *this;
  }

  // ---------------------------------------------------------------------------
  IncrementBase & IncrementBase::operator*=(const double & scalar)
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      auto view      = atlas::array::make_view<double, 2>(fields_[iv]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          view(n, l) *= scalar;
    }
    return *this;
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::axpy(const double & w, const IncrementBase & rhs)
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      auto       dst = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src = atlas::array::make_view<double, 2>(rhs.fields_[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) += w * src(n, l);
    }
  }

  // ---------------------------------------------------------------------------
  double IncrementBase::dot_product_with(const IncrementBase & other) const
  {
    const auto ghost = atlas::array::make_view<int32_t, 1>(
                           functionSpace_.ghost());
    double result = 0.0;
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      const auto v1      = atlas::array::make_view<double, 2>(fields_[name]);
      const auto v2      = atlas::array::make_view<double, 2>(
                               other.fields_[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
      {
        if (ghost(n)) continue;
        for (int l = 0; l < nlev; ++l)
          result += v1(n, l) * v2(n, l);
      }
    }
    comm_.allReduceInPlace(result, eckit::mpi::sum());
    return result;
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::schur_product_with(const IncrementBase & other)
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      auto       dst = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src = atlas::array::make_view<double, 2>(
                           other.fields_[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) *= src(n, l);
    }
  }

  // ---------------------------------------------------------------------------
  double IncrementBase::norm() const
  {
    const auto ghost = atlas::array::make_view<int32_t, 1>(
                           functionSpace_.ghost());
    double sumSq  = 0.0;
    size_t nOwned = 0;
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const auto view = atlas::array::make_view<double, 2>(fields_[iv]);
      const int npts  = fields_[iv].shape(0);
      const int nlev  = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
      {
        if (ghost(n)) continue;
        for (int l = 0; l < nlev; ++l)
        {
          sumSq += view(n, l) * view(n, l);
          ++nOwned;
        }
      }
    }
    comm_.allReduceInPlace(sumSq,  eckit::mpi::sum());
    comm_.allReduceInPlace(nOwned, eckit::mpi::sum());
    return nOwned > 0
        ? std::sqrt(sumSq / static_cast<double>(nOwned)) : 0.0;
  }

  // ---------------------------------------------------------------------------
  size_t IncrementBase::serialSize() const
  {
    const auto ghost = atlas::array::make_view<int32_t, 1>(
                           functionSpace_.ghost());
    size_t sz = time_.serialSize();
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        if (!ghost(n)) sz += static_cast<size_t>(nlev);
    }
    return sz;
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::serialize(std::vector<double> & vect) const
  {
    time_.serialize(vect);
    const auto ghost = atlas::array::make_view<int32_t, 1>(
                           functionSpace_.ghost());
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const auto view = atlas::array::make_view<double, 2>(fields_[iv]);
      const int npts  = fields_[iv].shape(0);
      const int nlev  = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
      {
        if (ghost(n)) continue;
        for (int l = 0; l < nlev; ++l)
          vect.push_back(view(n, l));
      }
    }
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::deserialize(const std::vector<double> & vect,
                                  size_t & current)
  {
    time_.deserialize(vect, current);
    const auto ghost = atlas::array::make_view<int32_t, 1>(
                           functionSpace_.ghost());
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      auto view      = atlas::array::make_view<double, 2>(fields_[iv]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
      {
        if (ghost(n)) continue;
        for (int l = 0; l < nlev; ++l)
          view(n, l) = vect[current++];
      }
    }
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::toFieldSet(atlas::FieldSet & fset) const
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
      fset.add(fields_[iv]);
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::fromFieldSet(const atlas::FieldSet & fset)
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      if (!fset.has(name)) continue;
      auto       dst = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src = atlas::array::make_view<double, 2>(fset[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      const int srcLev = fset[name].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) = (l < srcLev) ? src(n, l) : 0.0;
    }
  }

  // ---------------------------------------------------------------------------
  void IncrementBase::print(std::ostream & os) const
  {
    const auto ghost = atlas::array::make_view<int32_t, 1>(
                           functionSpace_.ghost());
    os << "Increment: time=" << time_ << "  nvars=" << vars_.size()
       << "  nlevels=" << numLevels_ << std::endl;
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      const auto view = atlas::array::make_view<double, 2>(fields_[iv]);
      const int npts  = fields_[iv].shape(0);
      const int nlev  = fields_[iv].shape(1);
      double mn = 0.0, mx = 0.0, rms = 0.0;
      size_t cnt = 0;
      for (int n = 0; n < npts; ++n)
      {
        if (ghost(n)) continue;
        for (int l = 0; l < nlev; ++l)
        {
          const double v = view(n, l);
          if (cnt == 0) { mn = mx = v; }
          mn   = std::min(mn, v);
          mx   = std::max(mx, v);
          rms += v * v;
          ++cnt;
        }
      }
      if (cnt > 0) rms = std::sqrt(rms / static_cast<double>(cnt));
      os << "  " << name << ": min=" << mn << "  max=" << mx
         << "  rms=" << rms << std::endl;
    }
  }

}  // namespace ijedi

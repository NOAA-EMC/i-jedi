// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include <algorithm>
#include <cmath>
#include <string>

#include "atlas/array.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/State/base/StateBase.h"
#include "ijedi/State/mom6/StateMOM6.h"

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

  std::shared_ptr<StateBase> StateBase::create(const Geometry & geom,
                                               const oops::Variables & vars,
                                               const util::DateTime & time)
  {
    const std::string type = geomType(geom);
    if (type == "mom6")
      return std::make_shared<StateMOM6>(geom, vars, time);
    throw eckit::BadValue("Unsupported geometry type for State: " + type,
                          Here());
  }

  std::shared_ptr<StateBase> StateBase::create(
      const Geometry & geom,
      const eckit::Configuration & conf)
  {
    const std::string type = geomType(geom);
    if (type == "mom6")
      return std::make_shared<StateMOM6>(geom, conf);
    throw eckit::BadValue("Unsupported geometry type for State: " + type,
                          Here());
  }

  std::shared_ptr<StateBase> StateBase::create(const Geometry & geom,
                                               const StateBase & other)
  {
    const std::string type = geomType(geom);
    if (type == "mom6")
      return std::make_shared<StateMOM6>(geom, other);
    throw eckit::BadValue("Unsupported geometry type for State: " + type,
                          Here());
  }

  std::shared_ptr<StateBase> StateBase::create(const oops::Variables & vars,
                                               const StateBase & other)
  {
    if (const auto * s = dynamic_cast<const StateMOM6 *>(&other))
      return std::make_shared<StateMOM6>(vars, *s);
    throw eckit::BadValue(
        "Unsupported StateBase subtype in variable-change constructor",
        Here());
  }

  // ---------------------------------------------------------------------------
  // Private helper — allocate one zero-filled field per variable
  // ---------------------------------------------------------------------------
  void StateBase::allocateFields()
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
  // Protected init helpers — called from model-specific subclass constructors
  // ---------------------------------------------------------------------------

  void StateBase::initFromGeomVarsTime(const Geometry & geom,
                                       const oops::Variables & vars,
                                       const util::DateTime & time)
  {
    functionSpace_ = geom.functionSpace();
    numLevels_     = geom.numLevels();
    vars_          = vars;
    time_          = time;
    allocateFields();
  }

  // ---------------------------------------------------------------------------
  void StateBase::initFromConf(const Geometry & geom,
                               const eckit::Configuration &conf)
  {
    functionSpace_ = geom.functionSpace();
    numLevels_     = geom.numLevels();
    vars_          = oops::Variables(conf.getStringVector("variables"));
    time_          = util::DateTime(conf.getString("date"));
    allocateFields();
  }

  // ---------------------------------------------------------------------------
  // Change-geometry: re-allocate fields on the new geometry and copy source
  // data for matching variables (point-by-point up to the smaller dimension).
  // When both geometries have the same layout this is an exact copy.
  void StateBase::initFromGeomState(const Geometry & geom,
                                    const StateBase & other)
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
      auto       dst  = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src  = atlas::array::make_view<double, 2>(
                            other.fields_[name]);
      const int npts  = std::min(fields_[iv].shape(0),
                                 other.fields_[name].shape(0));
      const int nlev  = std::min(fields_[iv].shape(1),
                                 other.fields_[name].shape(1));
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) = src(n, l);
    }
  }

  // ---------------------------------------------------------------------------
  // Change-variables: allocate fields for the new variable list, copying data
  // for variables that exist in the source and zeroing the rest.
  void StateBase::initFromVarsState(const Geometry & geom,
                                    const oops::Variables & vars,
                                    const StateBase & other)
  {
    functionSpace_ = geom.functionSpace();
    numLevels_     = geom.numLevels();
    vars_          = vars;
    time_          = other.time_;
    allocateFields();

    for (int iv = 0; iv < static_cast<int>(vars_.size()); ++iv)
    {
      const std::string & name = vars_[iv].name();
      if (!other.fields_.has(name)) continue;
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
  // Copy data (vars, time, field values) from another state of the same
  // geometry.  Does not touch functionSpace_, numLevels_, or comm_.
  void StateBase::copyDataFrom(const StateBase & other)
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
  void StateBase::zero()
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      auto view = atlas::array::make_view<double, 2>(fields_[iv]);
      view.assign(0.0);
    }
  }

  // ---------------------------------------------------------------------------
  void StateBase::accumul(const double & w, const StateBase & other)
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
          dst(n, l) += w * src(n, l);
    }
  }

  // ---------------------------------------------------------------------------
  double StateBase::norm() const
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
  size_t StateBase::serialSize() const
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
  void StateBase::serialize(std::vector<double> & vect) const
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
  void StateBase::deserialize(const std::vector<double> & vect,
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
  void StateBase::toFieldSet(atlas::FieldSet & fset) const
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
      fset.add(fields_[iv]);
  }

  // ---------------------------------------------------------------------------
  void StateBase::fromFieldSet(const atlas::FieldSet & fset)
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      const std::string & name = vars_[iv].name();
      if (!fset.has(name)) continue;
      auto       dst = atlas::array::make_view<double, 2>(fields_[name]);
      const auto src = atlas::array::make_view<double, 2>(fset[name]);
      const int npts = fields_[iv].shape(0);
      const int nlev = fields_[iv].shape(1);
      for (int n = 0; n < npts; ++n)
        for (int l = 0; l < nlev; ++l)
          dst(n, l) = src(n, l);
    }
  }

  // ---------------------------------------------------------------------------
  void StateBase::print(std::ostream & os) const
  {
    const auto ghost = atlas::array::make_view<int32_t, 1>(
                           functionSpace_.ghost());

    os << "State: time=" << time_ << "  nvars=" << vars_.size()
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

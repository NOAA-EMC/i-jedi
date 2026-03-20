// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "ijedi/State/mom6/StateMOM6.h"

namespace ijedi
{

  StateMOM6::StateMOM6(const Geometry & geom, const oops::Variables & vars,
                       const util::DateTime & time)
    : StateBase(geom.getComm()), geom_(geom)
  { initFromGeomVarsTime(geom, vars, time); }

  StateMOM6::StateMOM6(const Geometry & geom,
                       const eckit::Configuration & conf)
    : StateBase(geom.getComm()), geom_(geom)
  {
    initFromConf(geom, conf);
    read(conf);
  }

  StateMOM6::StateMOM6(const Geometry & geom, const StateBase & other)
    : StateBase(geom.getComm()), geom_(geom)
  { initFromGeomState(geom, other); }

  StateMOM6::StateMOM6(const oops::Variables & vars,
                       const StateMOM6 & other)
    : StateBase(other.comm_), geom_(other.geom_)
  { initFromVarsState(other.geom_, vars, other); }

  void StateMOM6::read(const eckit::Configuration & conf)
  {
    if (conf.getBool("analytic init", false)) { fillAnalytic_(); return; }
    if (!conf.has("filepath")) { return; }  // no I/O: fields remain zero

    const std::string fmt = conf.getString("file format", "fieldset");
    if (fmt == "fieldset") {
      oops::Log::trace() << "StateMOM6::read (fieldset) starting"
                         << std::endl;
      std::vector<std::string> varNames;
      std::vector<size_t>      varSizes;
      for (int iv = 0; iv < static_cast<int>(vars_.size()); ++iv) {
        varNames.push_back(vars_[iv].name());
        varSizes.push_back(static_cast<size_t>(numLevels_));
      }
      atlas::FieldSet fset;
      util::readFieldSet(comm_, functionSpace_, varSizes, varNames,
                         conf, fset);
      fromFieldSet(fset);
      oops::Log::trace() << "StateMOM6::read (fieldset) done" << std::endl;
    } else if (fmt == "native") {
      ABORT("StateMOM6::read — native MOM6 restart I/O not yet implemented");
    } else {
      throw eckit::BadValue(
          "StateMOM6::read — unknown file format: " + fmt, Here());
    }
  }

  void StateMOM6::write(const eckit::Configuration & conf) const
  {
    const std::string fmt = conf.getString("file format", "fieldset");
    if (fmt == "fieldset") {
      oops::Log::trace() << "StateMOM6::write (fieldset) starting"
                         << std::endl;
      util::writeFieldSet(comm_, conf, fields_);
      oops::Log::trace() << "StateMOM6::write (fieldset) done" << std::endl;
    } else if (fmt == "native") {
      ABORT(
          "StateMOM6::write — native MOM6 restart I/O not yet implemented");
    } else {
      throw eckit::BadValue(
          "StateMOM6::write — unknown file format: " + fmt, Here());
    }
  }

  void StateMOM6::fillAnalytic_()
  {
    for (int iv = 0; iv < fields_.size(); ++iv)
    {
      auto view = atlas::array::make_view<double, 2>(fields_[iv]);
      view.assign(1.0);
    }
  }

}  // namespace ijedi

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

#include "ijedi/Increment/mom6/IncrementMOM6.h"

namespace ijedi
{

  IncrementMOM6::IncrementMOM6(const Geometry & geom,
                               const oops::Variables & vars,
                               const util::DateTime & time)
    : IncrementBase(geom.getComm()), geom_(geom)
  { initFromGeomVarsTime(geom, vars, time); }

  IncrementMOM6::IncrementMOM6(const Geometry & geom,
                               const IncrementBase & other)
    : IncrementBase(geom.getComm()), geom_(geom)
  { initFromGeomIncrement(geom, other); }

  IncrementMOM6::IncrementMOM6(const IncrementMOM6 & other, bool copy)
    : IncrementBase(other.comm_), geom_(other.geom_)
  { initCopy(other, copy); }

  void IncrementMOM6::read(const eckit::Configuration & conf)
  {
    if (!conf.has("filepath")) { return; }  // no I/O: fields remain zero

    const std::string fmt = conf.getString("file format", "fieldset");
    if (fmt == "fieldset") {
      oops::Log::trace() << "IncrementMOM6::read (fieldset) starting"
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
      oops::Log::trace() << "IncrementMOM6::read (fieldset) done" << std::endl;
    } else if (fmt == "native") {
      ABORT("IncrementMOM6::read — native MOM6 I/O not yet implemented");
    } else {
      throw eckit::BadValue(
          "IncrementMOM6::read — unknown file format: " + fmt, Here());
    }
  }

  void IncrementMOM6::write(const eckit::Configuration & conf) const
  {
    const std::string fmt = conf.getString("file format", "fieldset");
    if (fmt == "fieldset") {
      oops::Log::trace() << "IncrementMOM6::write (fieldset) starting"
                         << std::endl;
      util::writeFieldSet(comm_, conf, fields_);
      oops::Log::trace() << "IncrementMOM6::write (fieldset) done"
                         << std::endl;
    } else if (fmt == "native") {
      ABORT("IncrementMOM6::write — native MOM6 I/O not yet implemented");
    } else {
      throw eckit::BadValue(
          "IncrementMOM6::write — unknown file format: " + fmt, Here());
    }
  }

}  // namespace ijedi

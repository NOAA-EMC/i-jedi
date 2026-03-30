// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include "ijedi/LinearVariableChange/LinearVariableChange.h"

#include "oops/util/Logger.h"

namespace ijedi
{

  LinearVariableChange::LinearVariableChange(const Geometry & /*geom*/,
                                             const eckit::Configuration & /*conf*/)
  {
    oops::Log::trace() << "LinearVariableChange::LinearVariableChange" << std::endl;
  }

}  // namespace ijedi

// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include "ijedi/ErrorCovariance/ErrorCovariance.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace ijedi
{

  ErrorCovariance::ErrorCovariance(const Geometry & /*geom*/,
                                   const oops::Variables & /*vars*/,
                                   const eckit::Configuration & /*conf*/,
                                   const State & /*xb*/,
                                   const State & /*fg*/)
  {
    throw eckit::NotImplemented(
        "ijedi::ErrorCovariance: use SABER for background error covariance",
        Here());
  }

  void ErrorCovariance::multiply(const Increment &, Increment &) const
  { throw eckit::NotImplemented("ijedi::ErrorCovariance::multiply", Here()); }

  void ErrorCovariance::inverseMultiply(const Increment &, Increment &) const
  { throw eckit::NotImplemented("ijedi::ErrorCovariance::inverseMultiply", Here()); }

  void ErrorCovariance::randomize(Increment &) const
  { throw eckit::NotImplemented("ijedi::ErrorCovariance::randomize", Here()); }

  void ErrorCovariance::print(std::ostream & os) const
  { os << "ErrorCovariance: stub (use SABER)"; }

}  // namespace ijedi

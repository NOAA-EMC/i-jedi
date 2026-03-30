// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class Geometry;
  class Increment;
  class State;

  // ---------------------------------------------------------------------------
  // Stub covariance required by oops::Traits.
  // The real background error covariance is provided by SABER.
  // All operational methods abort at runtime.
  // ---------------------------------------------------------------------------
  class ErrorCovariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ErrorCovariance>
  {
   public:
    static const std::string classname() { return "ijedi::ErrorCovariance"; }

    ErrorCovariance(const Geometry &, const oops::Variables &,
                    const eckit::Configuration &,
                    const State &, const State &);
    ~ErrorCovariance() = default;

    void multiply(const Increment &, Increment &) const;
    void inverseMultiply(const Increment &, Increment &) const;
    void randomize(Increment &) const;

   private:
    void print(std::ostream &) const override;
  };

}  // namespace ijedi

/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "mist/base/Increment.h"
#include "oops/util/ObjectCounter.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace oops {
class Variables;
}  // namespace oops

namespace util {
class DateTime;
}  // namespace util

namespace ijedi {

class Geometry;

class Increment : public mist::base::Increment, private util::ObjectCounter<Increment> {
 public:
  static std::string classname() { return "ijedi::Increment"; }

  Increment(const Geometry &, const oops::Variables &, const util::DateTime &);
  Increment(const Geometry &, const Increment &, const bool ad = false);
  Increment(const Increment &, const bool copy = true);
  ~Increment();

  Increment & operator=(const Increment &);

  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

 private:
  void print(std::ostream & os) const override;

  const Geometry & geom_;
};

}  // namespace ijedi

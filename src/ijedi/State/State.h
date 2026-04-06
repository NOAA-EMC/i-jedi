/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "mist/base/State.h"
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

class State : public mist::base::State, private util::ObjectCounter<State> {
 public:
  static std::string classname() { return "ijedi::State"; }

  State(const Geometry &, const eckit::Configuration &);
  State(const Geometry &, const oops::Variables &, const util::DateTime &, bool initToZero = true);
  State(const Geometry &, const State &);
  State(const oops::Variables &, const State &);
  State(const State &);

  ~State();

  State & operator=(const State &);

  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;

 private:
  void analytic_init(const eckit::Configuration &);

  void print(std::ostream & os) const override;

  const Geometry & geom_;
};

}  // namespace ijedi

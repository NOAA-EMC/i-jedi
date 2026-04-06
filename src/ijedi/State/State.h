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

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Io/IoBase.h"

namespace eckit
{
  class Configuration;
} // namespace eckit

namespace oops
{
  class Variables;
} // namespace oops

namespace util
{
  class DateTime;
} // namespace util

namespace ijedi
{

  class Geometry;

  // -------------------------------------------------------------------------------------------------

  class AnalyticICParameters : public oops::Parameters
  {
    OOPS_CONCRETE_PARAMETERS(AnalyticICParameters, Parameters)
  public:
    // Analytic initial condition parameters
    oops::RequiredParameter<std::string> method{"method", this};
  };

  // -------------------------------------------------------------------------------------------------

  class StateParameters : public oops::Parameters
  {
    OOPS_CONCRETE_PARAMETERS(StateParameters, Parameters)
  public:
    // Analytic initial condition parameters
    oops::OptionalParameter<oops::Variables> stateVariables{"state variables", this};
    oops::OptionalParameter<AnalyticICParameters> analytic{"analytic init", this};
    oops::OptionalParameter<util::DateTime> datetime{"datetime", this};
    // Read parameters
    IoParametersWrapper ioParametersWrapper{this};
    oops::OptionalParameter<bool> setdatetime{"set datetime on read", this};
  };

  // -------------------------------------------------------------------------------------------------

  class StateWriteParameters : public oops::Parameters
  {
    OOPS_CONCRETE_PARAMETERS(StateWriteParameters, Parameters)
  public:
    oops::OptionalParameter<std::string> type{"type", this};
    oops::OptionalParameter<std::string> exp{"exp", this};
    oops::OptionalParameter<int> member{"member", this};
    oops::OptionalParameter<std::string> memberPattern{"member pattern", this};
    oops::OptionalParameter<util::DateTime> date{"date", this};
    oops::OptionalParameter<int> iteration{"iteration", this};
    oops::OptionalParameter<std::string> prefix{"prefix", this};
    oops::Parameter<bool> dateCols{"date colons", true, this};
    IoParametersWrapper ioParametersWrapper{this};
    // Additional formats to output
    oops::OptionalParameter<std::vector<IoParametersWrapper>>
        additionalIo{"additional output formats", this};
  };

  // -------------------------------------------------------------------------------------------------

  class State : public mist::base::State, private util::ObjectCounter<State>
  {
  public:
    static std::string classname() { return "jedimpas::State"; }

    State(const Geometry &, const eckit::Configuration &);
    State(const Geometry &, const oops::Variables &, const util::DateTime &, bool initToZero = true);
    State(const Geometry &, const State &);
    State(const oops::Variables &, const State &);
    State(const State &);

    ~State();

    State &operator=(const State &);

    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &) const;

  private:
    void analytic_init(const eckit::Configuration &);

    void print(std::ostream &os) const override;

    const Geometry &geom_;
  };

} // namespace ijedi

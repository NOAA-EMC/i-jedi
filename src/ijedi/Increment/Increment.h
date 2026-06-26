#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "mist/base/Increment.h"
#include "oops/util/ObjectCounter.h"

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Io/IoBase.h"

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

  // -----------------------------------------------------------------------------------------------

  class DiracParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(DiracParameters, Parameters)
   public:
    // lon/lat dirac specification (grid-agnostic: nearest owned node on any
    // functionspace is found via a KD-tree). Levels are 1-based.
    oops::RequiredParameter<std::vector<double>> lon{"lon", this};
    oops::RequiredParameter<std::vector<double>> lat{"lat", this};
    oops::RequiredParameter<std::vector<int>> level{"level", this};
    oops::RequiredParameter<std::vector<std::string>> variable{"variable", this};
  };

  // -----------------------------------------------------------------------------------------------

  class IncrementParameters : public oops::Parameters
  {
    OOPS_CONCRETE_PARAMETERS(IncrementParameters, Parameters)
   public:
    // Io parameters wrapper for polymorphic IO parameters (nested under "io" key)
    oops::OptionalParameter<IoParametersWrapper> io{"io", this};
  };

  // -----------------------------------------------------------------------------------------------

  class IncrementWriteParameters : public oops::Parameters
  {
    OOPS_CONCRETE_PARAMETERS(IncrementWriteParameters, Parameters)
   public:
    // Io parameters for writing (nested under "io" key)
    oops::OptionalParameter<IoParametersWrapper> io{"io", this};
  };

  // -----------------------------------------------------------------------------------------------

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
    void dirac(const eckit::Configuration &);

   private:
    void print(std::ostream & os) const override;

    const Geometry & geom_;
  };

// -----------------------------------------------------------------------------------------------

}  // namespace ijedi

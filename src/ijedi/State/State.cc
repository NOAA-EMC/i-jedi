/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/State/State.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>

#include "atlas/field.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Utilities/PrintHelper.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/for_each.h"

#include "ijedi/Io/IoBase.h"

namespace ijedi
{

  State::State(const Geometry &geom, const eckit::Configuration &config)
      : mist::base::State(geom, oops::Variables(config, "state variables"),
                          util::DateTime(config.getString("date")), false),
        geom_(geom)
  {
    // If config has 'analytic init' then call analytic_init, else if config has 'io' then call read
    if (config.has("analytic init")) {
      analytic_init(config);
    } else if (config.has("io")) {
      read(config);
    } else {
      throw eckit::BadParameter("ijedi::State: config must have 'io' or 'analytic init'",
                                Here());
    }
  }

  State::State(const Geometry &geom, const oops::Variables &vars, const util::DateTime &time,
               bool initToZero)
      : mist::base::State(geom, vars, time, initToZero), geom_(geom) {}

  State::State(const Geometry &geom, const State &other)
      : mist::base::State(geom, other), geom_(geom) {}

  State::State(const oops::Variables &vars, const State &other)
      : mist::base::State(vars, other), geom_(other.geom_) {}

  State::State(const State &other) : mist::base::State(other), geom_(other.geom_) {}

  State::~State() = default;

  State &State::operator=(const State &rhs)
  {
    mist::base::State::operator=(rhs);
    return *this;
  }

  void State::read(const eckit::Configuration &config)
  {
    oops::Log::trace() << "ijedi::State::read starting" << std::endl;

    // Create a Parameters object
    StateParameters params;
    params.deserialize(config);

    // Check that there are IO parameters
    if (params.io.value() == boost::none ||
        params.io.value()->ioParameters.value() == nullptr)
    {
      throw eckit::BadParameter("ijedi::State::read: No IO parameters provided", Here());
    }

    // Get the polymorphic IO parameters
    const IoParametersBase &ioParams = *params.io.value()->ioParameters.value();

    // Create the IO object to use
    // ---------------------------
    std::unique_ptr<IoBase> io(IoFactory::create(geom_, ioParams));

    // Call read method of child
    // -------------------------
    io->readBase(this->fieldSet());

    oops::Log::trace() << "ijedi::State::read done" << std::endl;
  }

  void State::write(const eckit::Configuration &config) const
  {
    oops::Log::trace() << "ijedi::State::write starting" << std::endl;

    // Create a Parameters object
    StateWriteParameters params;
    params.deserialize(config);

    // Check that there are IO parameters
    if (params.io.value() == boost::none ||
        params.io.value()->ioParameters.value() == nullptr)
    {
      throw eckit::BadParameter("ijedi::State::write: No IO parameters provided", Here());
    }

    // Get the polymorphic IO parameters
    const IoParametersBase &ioParams = *params.io.value()->ioParameters.value();

    // Create the IO object to use
    // ---------------------------
    std::unique_ptr<IoBase> io(IoFactory::create(geom_, ioParams));

    // Call write method of child
    // --------------------------
    io->writeBase(this->fieldSet());

    oops::Log::trace() << "ijedi::State::write done" << std::endl;
  }

  void State::analytic_init(const eckit::Configuration &config)
  {
    oops::Log::trace() << "ijedi::State::analytic_init starting" << std::endl;
    // TODO(someone): implement analytic init
    oops::Log::trace() << "ijedi::State::analytic_init done" << std::endl;
  }

  void State::print(std::ostream &os) const
  {
    os << std::endl
       << "  Valid time: " << this->validTime() << std::endl
       << ", nFields = " << this->variables().size();

    const auto &comm = geom_.comm();
    const auto &fs = this->fieldSet();
    for (const auto &var : this->variables())
    {
      const atlas::Field &field = fs.field(var.name());
      const auto[globalMin, globalMax, rms] = fieldMinMaxRMS(comm, field);
      os << std::endl
         << var.name() << " : " << std::scientific << std::setprecision(16) << "Min=" << globalMin
         << ", Max=" << globalMax << ", RMS=" << rms;
    }
  }

}  // namespace ijedi

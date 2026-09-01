/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/python/PyModelBridge.h"

#include <memory>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "ijedi/Model/python/StubPyModelBridge.h"
#include "ijedi/Model/python/SubprocessPyModelBridge.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

std::unique_ptr<PyModelBridge> PyModelBridge::create(const eckit::Configuration & config,
                                                     const eckit::mpi::Comm & comm) {
  const std::string backend = config.getString("backend");

  if (backend == "stub") {
    return std::make_unique<StubPyModelBridge>(config);
  }

  if (backend == "subprocess") {
    return std::make_unique<SubprocessPyModelBridge>(config, comm);
  }

  if (backend == "embedded") {
    // Neither transport is built yet. Which of them to build is still open, and the answer
    // now depends on where JEDI's own state lives: with the state on the host, a separate
    // process costs one extra copy of a few tens of MB per oops step and buys clean memory
    // separation; with the state on a GPU, an in-process backend can pass a device pointer
    // through DLPack and copy nothing at all.
    throw eckit::NotImplemented("PyModelBridge: the 'embedded' backend is not built yet. "
                                "Use 'subprocess', which additionally keeps the model's "
                                "python environment separate from the one JEDI is built "
                                "against - they are not jointly satisfiable.", Here());
  }

  throw eckit::BadValue("PyModelBridge: unknown backend '" + backend +
                        "'; expected 'stub', 'embedded' or 'subprocess'", Here());
}

// -------------------------------------------------------------------------------------------------

void PyModelBridge::setTrajectory(const Fields &, const util::DateTime &) {
  throw eckit::NotImplemented("PyModelBridge: this backend declares no linear model, so it "
                              "has no trajectory to linearise about", Here());
}

// -------------------------------------------------------------------------------------------------

void PyModelBridge::advanceTL(const Fields &, Fields &, int) {
  throw eckit::NotImplemented("PyModelBridge: this backend provides no tangent linear model",
                              Here());
}

// -------------------------------------------------------------------------------------------------

void PyModelBridge::advanceAD(const Fields &, Fields &, int) {
  throw eckit::NotImplemented("PyModelBridge: this backend provides no adjoint model", Here());
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

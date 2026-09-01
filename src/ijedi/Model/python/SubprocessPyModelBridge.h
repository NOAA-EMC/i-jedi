/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "ijedi/Model/python/PyModelBridge.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace ijedi {

// -------------------------------------------------------------------------------------------------

/// \brief Drive a Python model in a separate process.
///
/// The model runs in its own interpreter, reached over one end of a socketpair. That
/// separation is not incidental: the Python environment these models need cannot be merged
/// with the one JEDI is built against. NeuralGCM requires numpy 2.x while spack-stack
/// provides 1.26, and running a JAX model and a PyTorch model in the same job would mean
/// resolving two dependency trees that are not jointly satisfiable. In separate processes the
/// question does not arise, and each model gets its own virtual environment.
///
/// A socket rather than the child's standard streams, because the frameworks write to stdout
/// and stderr freely and would corrupt a protocol sharing them. The child inherits stdout and
/// stderr from this process, so its chatter and any traceback land in the JEDI log.
///
/// Only the root task starts a server. Spawning one per MPI task would import the framework
/// and load the checkpoint once per rank - tens of seconds and gigabytes each, repeated for
/// nothing, since only the root marshals data. The other tasks receive the adapter's
/// declaration by broadcast so they can still validate their geometry against it, and their
/// bridge operations are no-ops.
///
/// The server is started once and lives for the whole run. Restarting it per step would pay
/// framework import and XLA compilation every time, which dominates everything else here -
/// the data movement this design is sometimes criticised for is a few tens of MB per output
/// step, against tens of seconds of compilation avoided.
class SubprocessPyModelBridge : public PyModelBridge {
 public:
  static const std::string classname() { return "ijedi::SubprocessPyModelBridge"; }

  SubprocessPyModelBridge(const eckit::Configuration & config, const eckit::mpi::Comm & comm);
  ~SubprocessPyModelBridge() override;

  const Declaration & declaration() const override { return declaration_; }

  void encode(const Fields & inputs, const util::DateTime & validTime) override;
  void advance(int steps) override;
  void decode(Fields & outputs) const override;
  void reset() override;

  void setTrajectory(const Fields & inputs, const util::DateTime & validTime) override;
  void advanceTL(const Fields & dxIn, Fields & dxOut, int steps) override;
  void advanceAD(const Fields & dxOut, Fields & dxIn, int steps) override;

 private:
  /// Send a request and return the response header, throwing with the Python traceback
  /// attached if the server reported a failure.
  std::string call(const std::string & requestJson, const Fields & send,
                   Fields * receive = nullptr) const;

  void writeAll(const void * data, size_t bytes) const;
  void readAll(void * data, size_t bytes) const;

  void start(const eckit::Configuration & config);
  void stop();

  /// Share the adapter's declaration, and any failure to obtain it, with the other tasks.
  /// A failure has to be broadcast too: if only the root threw, everyone else would wait
  /// forever in the broadcast that follows.
  std::string shareDeclaration(const eckit::mpi::Comm & comm, const std::string & json,
                               bool rootSucceeded, const std::string & rootError);

  bool isRoot_ = true;

  std::string python_;
  std::string script_;
  std::string adapter_;

  int socket_ = -1;
  int childPid_ = -1;

  Declaration declaration_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

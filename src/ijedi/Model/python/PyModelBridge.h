/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace eckit {
class Configuration;
}  // namespace eckit

namespace ijedi {

// -------------------------------------------------------------------------------------------------

/// \brief The seam between i-jedi (C++, MPI) and a model written in Python.
///
/// Nothing here is specific to one model or to one array framework. NeuralGCM is JAX, AIGFS
/// and a PyTorch physical model are others, and the difference between them belongs on the
/// Python side of this boundary, in an adapter (see ijedi_model_adapter.py). Adding a model
/// should be a Python-only change.
///
/// The interface is stateful rather than a set of pure functions, because the encoded model
/// state is a large framework-native object - a JAX pytree, a torch tensor tree - that must
/// stay on the Python side and stay on the device. It is never serialised across this
/// boundary. Only the dense fields at the ends of a leg cross, once per oops step rather than
/// once per the model's own (much shorter) internal timestep.
///
/// Fields crossing the boundary are global rather than MPI-partitioned, carry the adapter's
/// own variable names and units, and are laid out level-major over the function space's global
/// node ordering: value (k, n) sits at index k * nGlobalNodes + n. For a structured grid the
/// node ordering is row-major in (latitude, longitude), so the Python side can reshape the
/// buffer straight to (level, latitude, longitude) with no further work.
///
/// MPI. Every task constructs a bridge and every task gets the declaration, so that
/// configuration errors are found on all tasks at once rather than deadlocking on one. Only
/// the root task's bridge is asked to encode, advance or decode. A backend that is expensive
/// to bring up - importing a framework, loading weights - is free to do that lazily and only
/// where it is needed, as long as declaration() works everywhere.
class PyModelBridge {
 public:
  typedef std::map<std::string, std::vector<double>> Fields;

  /// \brief What the adapter says it needs and provides. This is how model knowledge reaches
  ///        the C++ side without the C++ side knowing about any particular model: the adapter
  ///        declares, and ModelPython checks the declaration against the geometry.
  struct Declaration {
    /// Human-readable identification, for logs.
    std::string description;
    /// The model's own internal timestep. The oops step must be a whole multiple of it.
    util::Duration timestep;
    /// Pressure levels the model is discretised on, in Pa, in the adapter's own order.
    /// Empty means the adapter imposes no vertical requirement.
    std::vector<double> levelsPa;
    /// Prognostic variables, keyed by the adapter's name, valued by the JEDI long name the
    /// adapter suggests it corresponds to. An empty value means the adapter has no suggestion
    /// and the yaml must supply the mapping.
    std::map<std::string, std::string> inputVariables;
    /// Variables the model needs supplied but does not predict, same keying.
    std::map<std::string, std::string> forcingVariables;
    /// Whether the adapter provides tangent-linear and adjoint operators.
    bool supportsLinear = false;
  };

  /// \brief Build the backend named by the "backend" key of \p config.
  static std::unique_ptr<PyModelBridge> create(const eckit::Configuration & config,
                                               const eckit::mpi::Comm & comm);

  virtual ~PyModelBridge() = default;

  virtual const Declaration & declaration() const = 0;

  // --- nonlinear ---------------------------------------------------------------------------

  /// \brief Build the internal state from dense fields. Replaces any state already held.
  virtual void encode(const Fields & inputs, const util::DateTime & validTime) = 0;

  /// \brief Advance the internal state by \p steps of the model's internal timestep.
  virtual void advance(int steps) = 0;

  /// \brief Read dense fields back out of the internal state.
  virtual void decode(Fields & outputs) const = 0;

  /// \brief Drop the internal state, keeping the model itself loaded. The model outlives a
  ///        forecast leg; the encoded state does not.
  virtual void reset() = 0;

  // --- tangent linear and adjoint ----------------------------------------------------------
  //
  // Reserved rather than implemented. Both JAX (jax.jvp / jax.vjp) and PyTorch
  // (torch.func.jvp / torch.func.vjp) derive these automatically, so for a differentiable
  // Python model they are a wrapping exercise rather than the multi-year effort they are for a
  // hand-coded model. Two things to respect when they are filled in: the operator being
  // linearised is the whole composition decode-advance-encode, not advance alone, because the
  // encoder and decoder are themselves learned; and for a stochastic model the random draw
  // must be frozen between the trajectory and the linear runs, or the adjoint will be
  // inconsistent while still passing a dot-product test.
  //
  // The signatures mirror oops: the perturbation is handed in and handed back each step,
  // because oops holds the Increment between calls to stepTL and stepAD.

  /// \brief Establish the nonlinear trajectory the linear operators linearise about.
  virtual void setTrajectory(const Fields & inputs, const util::DateTime & validTime);

  /// \brief Tangent linear: advance a perturbation \p steps forward about the trajectory.
  virtual void advanceTL(const Fields & dxIn, Fields & dxOut, int steps);

  /// \brief Adjoint of advanceTL: propagate a sensitivity \p steps backward.
  virtual void advanceAD(const Fields & dxOut, Fields & dxIn, int steps);
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

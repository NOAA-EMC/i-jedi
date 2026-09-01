"""The contract a Python forecast model implements to be driven by i-jedi.

i-jedi has one C++ class for Python models, ``ijedi::ModelPython``, and one transport,
``ijedi::PyModelBridge``. Everything that distinguishes one model from another lives here, on
the Python side, in an adapter. Adding a model - NeuralGCM, AIGFS, a PyTorch physical model -
is a Python-only change: subclass ``IjediModelAdapter``, register it, and name it in the yaml
under ``adapter:``. No C++ and no rebuild.

Framework independence
----------------------
Nothing in this contract assumes JAX or PyTorch. Arrays crossing the boundary are exchanged
through DLPack wherever the framework supports it (``jax.dlpack``, ``torch.utils.dlpack``),
which is what makes a zero-copy path possible once JEDI's own state lives on the GPU. Adapters
should therefore avoid forcing arrays through host numpy unless they have to.

Array layout
------------
Fields arrive and leave as one flat buffer per variable, laid out level-major over the
function space's global node ordering: value ``(k, n)`` at index ``k * n_nodes + n``. For a
structured grid the node ordering is row-major in (latitude, longitude), so::

    field.reshape(n_levels, n_lat, n_lon)

is the model's array with no transposition. Latitudes follow the geometry's own order, which
for the Gaussian grids these models use runs north to south. An adapter that needs a different
order flips it here, where the model's conventions are known - not in the C++.

Units are SI, and variables carry the adapter's own names. The JEDI-side names are the
adapter's business only insofar as it suggests them in its declaration.

Parallelism and lifetime
------------------------
The adapter is only ever asked to encode, advance or decode on one task. It should hold the
encoded state itself, in framework-native form, on the device: that state is large and must
never be serialised back across the boundary. ``encode`` is called once per forecast leg,
``advance``/``decode`` once per oops output step, ``reset`` at the end of the leg. The model
itself - weights, compiled graph - outlives the leg and should be loaded once in ``__init__``,
because compilation and checkpoint loading dominate startup for both JAX and PyTorch.

Transport independence
----------------------
Every adapter is expected to run under every transport. i-jedi chooses one at runtime with the
yaml's ``backend:`` key - in process through an embedded interpreter, or out of process
through a model server - and an adapter should not be able to tell which it got. The two
serve genuinely different deployments: in process is the only way to pass a device pointer
without a copy once JEDI's own state lives on a GPU, while out of process is the only way to
run models whose dependency trees cannot be resolved together, which is the normal case once
a JAX model and a PyTorch model are both in play.

What that costs the adapter is that it must not assume it owns the process. Under the
embedded transport it is a guest inside a JEDI MPI rank, possibly one of several on a node.
So:

* Do not mutate global framework state - ``jax.config.update``, ``torch.set_default_device``,
  thread-count and device-visibility environment variables. In a dedicated server process
  that is harmless; inside an MPI rank it reconfigures the host application underneath JEDI.
  Set what you need at construction, scoped to your own objects.
* Do not let the framework pre-allocate the whole GPU. JAX's default of claiming most of
  device memory on first use will fight both JEDI's own allocations and any sibling rank.
* Do not spawn threads or processes on the assumption that the machine is idle.
* Do not rely on interpreter shutdown to release anything. ``reset`` ends a forecast leg and
  the adapter may be reused many times in one run; free what you hold there.
* Do not read stdin or print to stdout. Under the out-of-process transport those may be the
  control channel. Use ``logging``.
* Do not force arrays through host numpy when the framework speaks DLPack - that discards the
  zero-copy path the embedded transport exists to provide.

In return the adapter can rely on: being constructed once and living for the whole run; being
called on one task only; and receiving ``encode`` before any ``advance`` or ``decode``, and
``reset`` before any subsequent ``encode``.

The practical test of all this is that the same yaml with ``backend:`` flipped must produce
identical output.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Dict, Mapping, MutableMapping

import numpy as np

Fields = MutableMapping[str, np.ndarray]


@dataclass
class Declaration:
    """What the adapter needs and provides.

    i-jedi checks this against the configured geometry before the first forecast, so a
    mismatched grid, level set or variable list is a startup error rather than a silently
    wrong forecast. Be accurate here; it is the only thing standing between a transposed
    field and a plausible-looking result.
    """

    #: Human-readable identification, for logs.
    description: str

    #: The model's own internal timestep. The oops output step must be a whole multiple.
    timestep: timedelta

    #: Prognostic variables: the adapter's name -> the JEDI long name it corresponds to.
    #: Leave the value empty if there is no sensible suggestion and the yaml must decide.
    input_variables: Dict[str, str]

    #: Variables the model needs supplied but does not predict, keyed the same way.
    forcing_variables: Dict[str, str] = field(default_factory=dict)

    #: Pressure levels the model is discretised on, in Pa, in the order the model uses them.
    #: Empty means the adapter imposes no vertical requirement.
    levels_pa: list = field(default_factory=list)

    #: Whether advance_tl and advance_ad are implemented.
    supports_linear: bool = False


class IjediModelAdapter(abc.ABC):
    """Base class for a Python model driven by i-jedi.

    The constructor receives the ``adapter config`` block from the yaml verbatim, so an adapter
    defines its own configuration - checkpoint paths, random seeds, forcing files - without
    anything in C++ needing to know about them. Load the model here.

    Read the transport-independence rules in the module docstring before implementing this:
    the adapter runs both inside a JEDI MPI rank and alone in a server process, and must
    behave identically in both.
    """

    def __init__(self, config: Mapping):
        self.config = config

    # -- declaration ------------------------------------------------------------------------

    @property
    @abc.abstractmethod
    def declaration(self) -> Declaration:
        """Describe the model. Called on every task, so it must be cheap and side-effect free."""

    # -- nonlinear --------------------------------------------------------------------------

    @abc.abstractmethod
    def encode(self, inputs: Fields, valid_time: datetime) -> None:
        """Build the internal state from dense fields, replacing any state already held."""

    @abc.abstractmethod
    def advance(self, steps: int) -> None:
        """Advance the internal state by ``steps`` internal timesteps.

        Prefer the framework's batched rollout - ``model.unroll`` under JAX, a compiled loop
        under PyTorch - over calling a single-step function in a Python loop. The whole point
        of taking ``steps`` rather than being called repeatedly is that the loop stays inside
        the compiled region.
        """

    @abc.abstractmethod
    def decode(self, outputs: Fields) -> None:
        """Write dense fields out of the internal state into ``outputs``.

        Only variables actually returned need to be set; i-jedi leaves anything absent
        untouched in the state rather than zeroing it.
        """

    @abc.abstractmethod
    def reset(self) -> None:
        """Drop the internal state. The model itself stays loaded for the next leg."""

    # -- tangent linear and adjoint ---------------------------------------------------------
    #
    # Optional, and the reason a differentiable Python model is interesting for data
    # assimilation in the first place: both JAX (jax.jvp / jax.vjp) and PyTorch
    # (torch.func.jvp / torch.func.vjp) derive these automatically, so what is a multi-year
    # hand-coding effort for a conventional model is a wrapping exercise here.
    #
    # Two traps when implementing them. The operator being linearised is the whole composition
    # decode-advance-encode, not advance alone, because the encoder and decoder are themselves
    # learned. And for a stochastic model the random draw must be frozen between the trajectory
    # and the linear runs - linearising about a different draw gives a wrong adjoint that will
    # still pass a dot-product test, because jvp and vjp are exact transposes of each other
    # whatever they are transposing.

    def set_trajectory(self, inputs: Fields, valid_time: datetime) -> None:
        """Establish the nonlinear trajectory the linear operators linearise about."""
        raise NotImplementedError(
            f"{type(self).__name__} declares no linear model")

    def advance_tl(self, dx_in: Fields, dx_out: Fields, steps: int) -> None:
        """Tangent linear: advance a perturbation forward about the trajectory."""
        raise NotImplementedError(
            f"{type(self).__name__} provides no tangent linear model")

    def advance_ad(self, dx_out: Fields, dx_in: Fields, steps: int) -> None:
        """Adjoint of advance_tl: propagate a sensitivity backward."""
        raise NotImplementedError(
            f"{type(self).__name__} provides no adjoint model")


# ---------------------------------------------------------------------------------------------
# Registry. The yaml's "adapter" key is looked up here.
# ---------------------------------------------------------------------------------------------

_ADAPTERS: Dict[str, type] = {}


def register(name: str):
    """Class decorator registering an adapter under the name the yaml will use."""

    def decorate(cls):
        if name in _ADAPTERS:
            raise ValueError(f"adapter '{name}' is already registered")
        _ADAPTERS[name] = cls
        return cls

    return decorate


def create(name: str, config: Mapping) -> IjediModelAdapter:
    """Instantiate the adapter registered under ``name``."""
    if name not in _ADAPTERS:
        raise KeyError(
            f"no adapter named '{name}'. Registered adapters: "
            f"{', '.join(sorted(_ADAPTERS)) or '(none)'}")
    return _ADAPTERS[name](config)


def registered() -> list:
    """Names of every registered adapter."""
    return sorted(_ADAPTERS)

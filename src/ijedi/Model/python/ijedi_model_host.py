"""Adapter dispatch, shared by every transport.

Both transports - the in-process one and the model server - drive an adapter through this
module rather than calling it directly. That is what makes "every adapter runs under every
transport" structural instead of aspirational: there is one place that knows how to look an
adapter up, build its declaration, marshal arrays in and out, and turn a Python failure into
something the C++ side can report. A second copy of that logic in a second transport would
drift, and adapters would start working under one and not the other.

Arrays arrive from C++ as flat float64 buffers, laid out level-major over the function
space's global node ordering (value (k, n) at index k * n_nodes + n). Reshaping them into
whatever the model wants is the adapter's job, because only the adapter knows the model's
conventions.
"""

from __future__ import annotations

import datetime
import importlib
import traceback
from typing import Dict, Mapping

import numpy as np

from ijedi_model_adapter import create, registered  # noqa: F401

# Adapter modules are imported for their registration side effect. Each must import its
# framework lazily, inside __init__, so that listing the available adapters does not require
# every model's dependencies to be installed.
_ADAPTER_MODULES = [
    "adapters.neuralgcm_adapter",
]


def _load_adapter_modules():
    problems = {}
    for name in _ADAPTER_MODULES:
        try:
            importlib.import_module(name)
        except Exception as exc:  # noqa: BLE001 - a broken adapter must not hide the others
            problems[name] = f"{type(exc).__name__}: {exc}"
    return problems


class Host:
    """Owns one adapter instance for the lifetime of a run."""

    def __init__(self, adapter_name: str, config: Mapping):
        problems = _load_adapter_modules()
        try:
            self.adapter = create(adapter_name, config)
        except KeyError as exc:
            detail = ""
            if problems:
                detail = ("\nSome adapter modules failed to import, which may be why:\n  " +
                          "\n  ".join(f"{k}: {v}" for k, v in problems.items()))
            raise KeyError(str(exc) + detail) from None

    # -- declaration ----------------------------------------------------------------------

    def declaration(self) -> Dict:
        """The adapter's declaration, flattened to JSON-safe types for the wire."""
        d = self.adapter.declaration
        return {
            "description": d.description,
            "timestep_seconds": int(d.timestep.total_seconds()),
            "levels_pa": [float(v) for v in d.levels_pa],
            "input_variables": dict(d.input_variables),
            "forcing_variables": dict(d.forcing_variables),
            "supports_linear": bool(d.supports_linear),
        }

    # -- nonlinear ------------------------------------------------------------------------

    def encode(self, fields: Dict[str, np.ndarray], valid_time: str) -> None:
        when = datetime.datetime.fromisoformat(valid_time.replace("Z", "+00:00"))
        self.adapter.encode(fields, when)

    def advance(self, steps: int) -> None:
        self.adapter.advance(int(steps))

    def decode(self) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        self.adapter.decode(out)
        return {k: np.ascontiguousarray(v, dtype=np.float64).ravel()
                for k, v in out.items()}

    def reset(self) -> None:
        self.adapter.reset()

    # -- dispatch -------------------------------------------------------------------------

    def call(self, op: str, header: Mapping, arrays: Dict[str, np.ndarray]) -> Dict:
        """Run one operation. Returns a response header; arrays go under 'arrays'."""
        if op == "declare":
            return {"status": "ok", "declaration": self.declaration()}
        if op == "encode":
            self.encode(arrays, header["valid_time"])
            return {"status": "ok"}
        if op == "advance":
            self.advance(header["steps"])
            return {"status": "ok"}
        if op == "decode":
            return {"status": "ok", "arrays": self.decode()}
        if op == "reset":
            self.reset()
            return {"status": "ok"}
        if op == "set_trajectory":
            when = datetime.datetime.fromisoformat(
                header["valid_time"].replace("Z", "+00:00"))
            self.adapter.set_trajectory(arrays, when)
            return {"status": "ok"}
        if op == "advance_tl":
            out = {}
            self.adapter.advance_tl(arrays, out, int(header["steps"]))
            return {"status": "ok", "arrays": {k: np.ascontiguousarray(v, dtype=np.float64)
                                               for k, v in out.items()}}
        if op == "advance_ad":
            out = {}
            self.adapter.advance_ad(arrays, out, int(header["steps"]))
            return {"status": "ok", "arrays": {k: np.ascontiguousarray(v, dtype=np.float64)
                                               for k, v in out.items()}}
        raise ValueError(f"unknown operation '{op}'")


def error_response(exc: BaseException) -> Dict:
    """Render a failure so the C++ side can report it usefully rather than just 'failed'."""
    return {
        "status": "error",
        "message": f"{type(exc).__name__}: {exc}",
        "traceback": traceback.format_exc(),
    }

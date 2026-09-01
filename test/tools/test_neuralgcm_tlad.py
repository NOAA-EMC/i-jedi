#!/usr/bin/env python3
"""Check the NeuralGCM adapter's tangent linear against its adjoint.

The test is the dot-product identity <M dx, w> == <dx, M^T w>. For a linear model derived by
automatic differentiation this should hold to rounding, because jax.jvp and jax.vjp produce
exact transposes of one another - so a failure here means the two are being taken about
different operators, not that the derivative is inaccurate.

What this does NOT establish is that the tangent linear is a good approximation to the
nonlinear model over an assimilation window. That is a separate, and for a neural model
genuinely open, question: measured convergence of the finite difference is sub-first-order,
and the relative error at a perturbation of 1e-2 of the field standard deviation is around
30%. It is deliberately not asserted here, because there is no threshold that would be both
meaningful and passing, and a tolerance loosened until it passes tests nothing.

Usage: test_neuralgcm_tlad.py <initial_condition.nc> [--steps 1] [--tolerance 1e-5]
"""

import argparse
import datetime
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..",
    "src", "ijedi", "Model", "python"))

UPPER = ["u_component_of_wind", "v_component_of_wind", "temperature", "geopotential",
         "specific_humidity", "specific_cloud_liquid_water_content",
         "specific_cloud_ice_water_content"]
SURFACE = ["sea_surface_temperature", "sea_ice_cover"]


def load_fields(path):
    """Read the initial condition into i-jedi's flat layout.

    The file has latitude ascending; atlas numbers latitudes from north to south, so the
    array is flipped before flattening. This mirrors what the C++ gather produces.
    """
    import netCDF4 as nc
    dataset = nc.Dataset(path)
    fields = {}
    for name in UPPER:
        a = np.asarray(dataset[name][0]).astype(np.float64)
        fields[name] = np.ascontiguousarray(a[:, ::-1, :]).ravel()
    for name in SURFACE:
        a = np.asarray(dataset[name][0]).astype(np.float64)
        fields[name] = np.ascontiguousarray(a[::-1, :]).ravel()
    return fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("initial_condition")
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--tolerance", type=float, default=1.0e-5)
    args = parser.parse_args()

    from ijedi_model_host import Host

    fields = load_fields(args.initial_condition)
    host = Host("neuralgcm", {"checkpoint": "demo", "rng seed": "0"})

    declaration = host.declaration()
    if not declaration["supports_linear"]:
        print("FAILED: the adapter declares no linear model")
        return 1
    print(f"adapter: {declaration['description']}")

    adapter = host.adapter
    adapter.set_trajectory(fields, datetime.datetime(1959, 1, 2))

    rng = np.random.default_rng(42)
    dx = {k: rng.standard_normal(fields[k].size) * float(np.std(fields[k])) * 1.0e-3
          for k in UPPER}

    dy = {}
    adapter.advance_tl(dx, dy, args.steps)

    w = {k: rng.standard_normal(dy[k].size) for k in dy}
    adj = {}
    adapter.advance_ad(w, adj, args.steps)

    lhs = sum(float(np.dot(dy[k], w[k])) for k in dy)
    rhs = sum(float(np.dot(dx[k], adj[k])) for k in adj if k in dx)
    relative = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0e-30)

    print(f"  <M dx, w>   = {lhs: .12e}")
    print(f"  <dx, M^T w> = {rhs: .12e}")
    print(f"  relative difference = {relative:.3e}  (tolerance {args.tolerance:.1e})")

    if relative > args.tolerance:
        print("FAILED: the adjoint is not the transpose of the tangent linear")
        return 1
    print("adjoint test passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

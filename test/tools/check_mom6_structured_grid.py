#!/usr/bin/env python3
"""
Compare mom6_structured_grid.nc against the SOCA reference gridspec.

Fields compared (tolerance accounts for CDL text roundtrip ~1e-13):
  lon    <- lon    (atol 1e-10)
  lat    <- lat    (atol 1e-10)
  dxT    <- dx     (rtol 1e-9)
  dyT    <- dy     (rtol 1e-9)
  areaT  <- area   (rtol, default 1e-9; see --rtol-area)
  mask2d <- mask2d (exact)
  lonu   <- lonu   (atol 1e-10)
  latu   <- latu   (atol 1e-10)
  lonv   <- lonv   (atol 1e-10)
  latv   <- latv   (atol 1e-10)

Note on areaT: SOCA's reference area is read directly from MOM6's runtime grid
(spherical-geometry initialization), while i-jedi sums the 4 supergrid sub-cells.
This produces a systematic ~0.8% difference on real grids; use --rtol-area=1e-2
for real-data (tier-2) tests.  Synthetic test grids match at machine precision.

Usage:
    python3 check_mom6_structured_grid.py <mom6_structured_grid.nc> \
        <soca_gridspec.nc> [--rtol-area=<tol>]
"""

import sys
import netCDF4 as nc
import numpy as np


def check(name, arr, ref, atol=0.0, rtol=0.0):
    diff = np.abs(arr - ref)
    if rtol > 0.0:
        scale = np.abs(ref) + 1e-30
        err = np.max(diff / scale)
        ok = err <= rtol
        label = f"max rel diff = {err:.3e}"
    else:
        err = np.max(diff)
        ok = err <= atol
        label = f"max abs diff = {err:.3e}"
    print(f"  {name:8s}  {label}  {'OK' if ok else 'FAIL'}")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = {k: v for k, v in (a[2:].split("=") for a in sys.argv[1:] if a.startswith("--"))}

    if len(args) != 2:
        print(
            f"Usage: {sys.argv[0]}"
            " <mom6_structured_grid.nc> <soca_gridspec.nc> [--rtol-area=<tol>]"
        )
        sys.exit(1)

    grid_file = args[0]
    soca_file = args[1]
    rtol_area = float(opts.get("rtol-area", "1e-9"))

    g = nc.Dataset(grid_file)
    s = nc.Dataset(soca_file)

    # SOCA gridspec has a leading time dimension of size 1
    soca = {
        "lon":    s.variables["lon"][0],
        "lat":    s.variables["lat"][0],
        "dx":     s.variables["dx"][0],
        "dy":     s.variables["dy"][0],
        "area":   s.variables["area"][0],
        "mask2d": s.variables["mask2d"][0],
        "lonu":   s.variables["lonu"][0],
        "latu":   s.variables["latu"][0],
        "lonv":   s.variables["lonv"][0],
        "latv":   s.variables["latv"][0],
    }

    def gv(name):
        return g.variables[name][:]

    passed = True
    print(f"Comparing\n  {grid_file}\nvs\n  {soca_file}")
    passed &= check("lon",    gv("lon"),    soca["lon"],    atol=1e-10)
    passed &= check("lat",    gv("lat"),    soca["lat"],    atol=1e-10)
    passed &= check("dxT",    gv("dxT"),    soca["dx"],     rtol=1e-9)
    passed &= check("dyT",    gv("dyT"),    soca["dy"],     rtol=1e-9)
    passed &= check("areaT",  gv("areaT"),  soca["area"],   rtol=rtol_area)
    passed &= check("mask2d", gv("mask2d"), soca["mask2d"], atol=0.0)
    passed &= check("lonu",   gv("lonu"),   soca["lonu"],   atol=1e-10)
    passed &= check("latu",   gv("latu"),   soca["latu"],   atol=1e-10)
    passed &= check("lonv",   gv("lonv"),   soca["lonv"],   atol=1e-10)
    passed &= check("latv",   gv("latv"),   soca["latv"],   atol=1e-10)

    g.close()
    s.close()

    if passed:
        print("All checks passed.")
        sys.exit(0)
    else:
        print("One or more checks FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    main()

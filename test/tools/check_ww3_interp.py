#!/usr/bin/env python3
"""Check the MOM6 <-> WW3 convertstate outputs against the analytic fields.

Usage: check_ww3_interp.py <data_generated dir>

  ww3/ww3_to_mom6.nc  hs, aice on the MOM6 grid: checked at wet MOM6 cells
  ww3/mom6_to_ww3.nc  ssh, ice on the WW3 nodes: checked at every node

Values not reached by the interpolation (missing / _FillValue) are counted,
not checked, and must stay below MAX_MISSING of the points: ~12% of the wet
MOM6 cells, along the coasts, are outside the coarse WW3 test mesh.
"""

import os
import sys

import netCDF4
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_mom6_testdata import ice_concentration  # noqa: E402
from gen_ww3_testdata import wave_height  # noqa: E402

MAX_MISSING = 0.15
# max |error|: the WW3 test mesh is 4 degrees, the MOM6 grid 5. Ice is a
# ramp of 0.1/degree with kinks (e.g. at 70N, where WW3 nodes sit): linear
# interpolation across a kink from MOM6 cells 2.5 degrees away errs by 0.25
TOL = {"hs": 0.01, "aice": 0.15, "ssh": 0.01, "ice": 0.3}


def check(name, value, exact):
    value = np.asarray(value, float)
    missing = ~np.isfinite(value) | (np.abs(value) > 1e10)
    err = np.abs(value - exact)[~missing].max()
    frac = missing.mean()
    ok = err <= TOL[name] and frac <= MAX_MISSING
    print(f"{name:5s}: max |error| {err:.4f} (tol {TOL[name]}), "
          f"missing {missing.sum()}/{missing.size} {'OK' if ok else 'FAIL'}")
    return ok


def main():
    d = sys.argv[1]
    ok = True

    with netCDF4.Dataset(os.path.join(d, "mom6", "soca_gridspec.72x35x25.nc")) as gs:
        lon, lat = gs["lon"][0].data, gs["lat"][0].data
        wet = gs["mask2d"][0].data > 0.5
    with netCDF4.Dataset(os.path.join(d, "ww3", "ww3_to_mom6.nc")) as ds:
        ds.set_auto_mask(False)
        ok &= check("hs", ds["hs"][:].squeeze()[wet], wave_height(lon, lat)[wet])
        ok &= check("aice", ds["aice"][:].squeeze()[wet], ice_concentration(lat)[wet])

    with netCDF4.Dataset(os.path.join(d, "ww3", "mom6_to_ww3.nc")) as ds:
        ds.set_auto_mask(False)
        lon, lat = ds["lon"][:], ds["lat"][:]
        ok &= check("ssh", ds["ssh"][:], 0.1 * np.sin(np.deg2rad(lat)))
        ok &= check("ice", ds["ice"][:], ice_concentration(lat))

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

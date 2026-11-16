#!/usr/bin/env python3
"""Check the sea-ice io of i-jedi on the 72x35 MOM6 test grid.

Usage: check_cice_io.py <data_generated/mom6 dir>

Expects, in that directory (gen_mom6_testdata.py and the ijedi_state_cice /
ijedi_state_mom6_cice_history tests):
  cice6.nc                         CICE history (aice_h, grid-cell mean hi_h, hs_h)
  cice6_restart.nc                 CICE restart template (aicen, vicen, vsnon, Tsfcn, iceumask)
  cice6_restart_write_test.nc      restart written by the cice io (update mode)
  mom6_cice_history_write_test.nc  state read with "thickness from volume", written by the
                                   mom6 io
  soca_gridspec.72x35x25.nc        latitudes

The oops state test multiplies the state by 3 before writing it.

Checks
  1. cice io update mode: the variables JEDI holds are 3 x the template everywhere, and
     iceumask, which JEDI does not hold, is the template's.
  2. thickness from volume: where there is ice, the thicknesses written are 3 x the
     per-ice-area values of the generator (ice 2 m north / 1 m south, snow 0.3 m / 0.1 m),
     i.e. hi_h / aice_h and hs_h / aice_h.
"""

import os
import sys

import netCDF4
import numpy as np

STATE_FACTOR = 3.0
TOL = 1.0e-6


def read(path, name):
    with netCDF4.Dataset(path) as ds:
        return np.squeeze(np.ma.filled(ds[name][:], np.nan)).astype(np.float64)


def check_restart(datadir):
    template = os.path.join(datadir, "cice6_restart.nc")
    written = os.path.join(datadir, "cice6_restart_write_test.nc")
    ok = True
    for name in ("aicen", "vicen", "vsnon", "Tsfcn"):
        err = np.nanmax(np.abs(read(written, name) - STATE_FACTOR * read(template, name)))
        print(f"  restart {name}: max |written - 3 x template| = {err:.3e}")
        ok = ok and err < TOL
    if not np.array_equal(read(written, "iceumask"), read(template, "iceumask")):
        print("  restart iceumask: changed by the write")
        ok = False
    else:
        print("  restart iceumask: unchanged")
    return ok


def check_history(datadir):
    history = os.path.join(datadir, "cice6.nc")
    written = os.path.join(datadir, "mom6_cice_history_write_test.nc")
    with netCDF4.Dataset(os.path.join(datadir, "soca_gridspec.72x35x25.nc")) as gs:
        lat = gs["lat"][0].data
    aice = read(history, "aice_h")
    ice = aice > 0.0
    if not ice.any():
        print("  history: no ice in the test data")
        return False
    ok = True
    for name, north, south in (("hi_h", 2.0, 1.0), ("hs_h", 0.3, 0.1)):
        expected = STATE_FACTOR * np.where(lat > 0.0, north, south)
        err = np.max(np.abs(read(written, name)[ice] - expected[ice]))
        print(f"  history {name} / aice_h: max error = {err:.3e} ({ice.sum()} ice points)")
        ok = ok and err < 1.0e-5
    return ok


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    datadir = sys.argv[1]
    ok = check_restart(datadir)
    ok = check_history(datadir) and ok
    print("PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

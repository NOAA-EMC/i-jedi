#!/usr/bin/env python3
"""
Generate test data for the doubly periodic Cartesian (atlas regional) grid.

The grid is a small patch centred on the equator, where lon/lat in degrees is
close to Cartesian x/y. The grid constants below must match the `grid:` block
of the tier1_integration_*-atlas_bifourier.yaml tests.

Outputs:
    1) background.nc: air_temperature on the grid, constant on each level, in the
       single-file StructuredColumns layout read by the "unstructured" io
    2) obs.nc: IODA airTemperature observations located exactly at grid points,
       including the domain corner and edges so that the analysis increment has
       to wrap around the periodic boundaries

Usage: python3 gen_cartesian_testdata.py <output_dir>
"""

import calendar
import datetime
import os
import sys

import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    print("netCDF4 python module is required", file=sys.stderr)
    sys.exit(1)

# Grid (must match the YAML `grid:` block)
NX = 64
NY = 48
NZ = 5
DLON = 0.09  # degrees (~10 km)
DLAT = 0.09  # degrees (~10 km)
WEST = -0.5 * (NX - 1) * DLON
SOUTH = -0.5 * (NY - 1) * DLAT

# Background: constant on each level, levels top down
BKG_PROFILE = [220.0, 235.0, 250.0, 265.0, 280.0]

# Observations at grid points (i, j), 0-based with j = 0 at the southern edge:
# south-west corner, middle of the west edge, middle of the south edge, centre
OBS_POINTS = [(0, 0), (0, NY // 2), (NX // 2, 0), (NX // 2, NY // 2)]
OBS_INNOVATION = 1.0
OBS_ERROR = 1.0
OBS_DATE = datetime.datetime(2020, 12, 15, 0, 0, 0)


def grid_lon(i):
    return WEST + i * DLON


def grid_lat(j):
    return SOUTH + j * DLAT


def generate_background(filepath):
    ds = nc.Dataset(filepath, "w", format="NETCDF4")
    ds.createDimension("nx", NX)
    ds.createDimension("ny", NY)
    ds.createDimension("nz_air_temperature", NZ)
    v = ds.createVariable("air_temperature", "f8", ("nz_air_temperature", "ny", "nx"))
    data = np.empty((NZ, NY, NX))
    for k in range(NZ):
        data[k, :, :] = BKG_PROFILE[k]
    v[:] = data
    ds.close()


def generate_obs(filepath):
    nlocs = len(OBS_POINTS)
    ds = nc.Dataset(filepath, "w", format="NETCDF4")
    ds.setncattr("_ioda_layout", "ObsGroup")
    ds.setncattr("_ioda_layout_version", np.int32(0))
    ds.createDimension("Location", nlocs)
    v = ds.createVariable("Location", "i4", ("Location",))
    v[:] = np.zeros(nlocs, dtype=np.int32)

    meta = ds.createGroup("MetaData")
    v = meta.createVariable("latitude", "f4", ("Location",))
    v.units = "degrees_north"
    v[:] = np.array([grid_lat(j) for (_, j) in OBS_POINTS], dtype=np.float32)
    v = meta.createVariable("longitude", "f4", ("Location",))
    v.units = "degrees_east"
    v[:] = np.array([grid_lon(i) for (i, _) in OBS_POINTS], dtype=np.float32)
    v = meta.createVariable("dateTime", "i8", ("Location",))
    v.units = "seconds since 1970-01-01T00:00:00Z"
    v[:] = np.full(nlocs, calendar.timegm(OBS_DATE.timetuple()), dtype=np.int64)

    # Identity operator uses the lowest model level
    obsval = ds.createGroup("ObsValue")
    v = obsval.createVariable("airTemperature", "f4", ("Location",))
    v.units = "K"
    v[:] = np.full(nlocs, BKG_PROFILE[-1] + OBS_INNOVATION, dtype=np.float32)

    obserr = ds.createGroup("ObsError")
    v = obserr.createVariable("airTemperature", "f4", ("Location",))
    v.units = "K"
    v[:] = np.full(nlocs, OBS_ERROR, dtype=np.float32)
    ds.close()


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output_dir>")
        sys.exit(1)

    outdir = sys.argv[1]
    print("Generating Cartesian grid test data in:", outdir)
    os.makedirs(outdir, exist_ok=True)

    bkg_path = os.path.join(outdir, "background.nc")
    generate_background(bkg_path)
    print(f"  {bkg_path}")

    obs_path = os.path.join(outdir, "obs.nc")
    generate_obs(obs_path)
    print(f"  {obs_path}")


if __name__ == "__main__":
    main()

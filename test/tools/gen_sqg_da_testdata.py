#!/usr/bin/env python3
"""
Generate a truth, observations and an ensemble for the SQG data assimilation tests.

The SQG domain is the doubly periodic 64x64 equatorial lon/lat patch of the SQG tests
(0.5 deg spacing, x -> lon 0..31.5E, y -> lat 15.75S..15.75N); the constants below must match
the `grid:` block of tier1_integration_3dvar-sqg.yaml and tier1_integration_letkf-sqg.yaml.

Outputs (in <data_dir>, which must contain sqg_ic.nc from gen_sqg_testdata.py):
    1) sqg_truth.nc: the initial condition plus a random perturbation on each boundary, drawn
       from the doubly periodic GC99 correlation used by the background error covariance (same
       support), with PERT_STD standard deviation; SQG file layout (pv in m/s)
    2) obs_sqg_bottom.nc and obs_sqg_lid.nc: IODA potentialTemperature observations of the
       truth plus random noise, at NOBS randomly selected grid points on the lower boundary
       and on the lid
    3) ens_memNNN.nc: NENS ensemble members, the initial condition plus random perturbations
       drawn like the truth perturbation and recentred on the initial condition; SQG layout

The random seeds are fixed so the output is reproducible.

Usage: python3 gen_sqg_da_testdata.py <data_dir>
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
N = 64
DLON = 0.5            # degrees
LON0 = 0.0
LAT0 = -15.75
DX = DLON * np.pi / 180.0 * 6371229.0   # grid spacing on the patch (m)

# Truth perturbation: GC99 support (m) and standard deviation (K), as in the B
PERT_SUPPORT = 400.0e3
PERT_STD = 1.0

# Observations on each boundary
NOBS = 200
OBS_ERROR = 0.5       # K
OBS_DATE = datetime.datetime(2015, 1, 1, 0, 0, 0)
SEED = 20261014

# Ensemble (perturbations drawn with their own seed, so the truth and obs do not depend on it)
NENS = 10
ENS_SEED = 20261015


def gc99(r):
    r = np.abs(r)
    out = np.zeros_like(r)
    m1 = r < 0.5
    m2 = (r >= 0.5) & (r < 1.0)
    x = r[m1]
    out[m1] = -8.0*x**5 + 8.0*x**4 + 5.0*x**3 - 20.0/3.0*x**2 + 1.0
    x = r[m2]
    out[m2] = 8.0/3.0*x**5 - 8.0*x**4 + 5.0*x**3 + 20.0/3.0*x**2 - 10.0*x + 4.0 - 1.0/(3.0*x)
    return out


def periodic_gc99_sqrt_spectrum():
    # Correlation with the origin on the periodic grid (nearest periodic image distances)
    i = np.minimum(np.arange(N), N - np.arange(N))
    dist = np.hypot(i[None, :] * DX, i[:, None] * DX)
    spectrum = np.fft.fft2(gc99(dist / PERT_SUPPORT)).real
    return np.sqrt(np.maximum(spectrum, 0.0))


def write_sqg(path, ic, theta, comment):
    scale = float(ic["pv"].getncattr("scale_to_K"))
    ds = nc.Dataset(path, "w", format="NETCDF4")
    for name in ("t", "z", "y", "x"):
        ds.createDimension(name, None if name == "t" else len(ic.dimensions[name]))
    v = ds.createVariable("pv", "f4", ("t", "z", "y", "x"))
    v.units = "m/s"
    v.scale_to_K = np.float32(scale)
    v[0] = (theta / scale).astype(np.float32)
    ds.setncatts({a: ic.getncattr(a) for a in ic.ncattrs() if a != "comment"})
    ds.comment = comment
    ds.close()


def write_obs(path, rng, truth_level):
    points = rng.choice(N * N, size=NOBS, replace=False)
    iy, ix = np.divmod(points, N)
    values = truth_level[iy, ix] + OBS_ERROR * rng.standard_normal(NOBS)

    ds = nc.Dataset(path, "w", format="NETCDF4")
    ds.setncattr("_ioda_layout", "ObsGroup")
    ds.setncattr("_ioda_layout_version", np.int32(0))
    ds.createDimension("Location", NOBS)
    v = ds.createVariable("Location", "i4", ("Location",))
    v[:] = np.zeros(NOBS, dtype=np.int32)

    meta = ds.createGroup("MetaData")
    v = meta.createVariable("latitude", "f4", ("Location",))
    v.units = "degrees_north"
    v[:] = (LAT0 + DLON * iy).astype(np.float32)
    v = meta.createVariable("longitude", "f4", ("Location",))
    v.units = "degrees_east"
    v[:] = (LON0 + DLON * ix).astype(np.float32)
    v = meta.createVariable("dateTime", "i8", ("Location",))
    v.units = "seconds since 1970-01-01T00:00:00Z"
    v[:] = np.full(NOBS, calendar.timegm(OBS_DATE.timetuple()), dtype=np.int64)

    obsval = ds.createGroup("ObsValue")
    v = obsval.createVariable("potentialTemperature", "f4", ("Location",))
    v.units = "K"
    v[:] = values.astype(np.float32)

    obserr = ds.createGroup("ObsError")
    v = obserr.createVariable("potentialTemperature", "f4", ("Location",))
    v.units = "K"
    v[:] = np.full(NOBS, OBS_ERROR, dtype=np.float32)
    ds.close()


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <data_dir>")
        sys.exit(1)

    datadir = sys.argv[1]
    print("Generating SQG 3DVar test data in:", datadir)

    # Background: initial condition, in K
    ic = nc.Dataset(os.path.join(datadir, "sqg_ic.nc"))
    scale = float(ic["pv"].getncattr("scale_to_K"))
    background = np.array(ic["pv"][0], dtype=np.float64) * scale   # (z, y, x)
    assert background.shape == (2, N, N)

    # Truth: background plus independent periodic GC99 perturbations on each boundary
    rng = np.random.default_rng(SEED)
    sqrt_spectrum = periodic_gc99_sqrt_spectrum()
    truth = background.copy()
    for k in range(2):
        noise = rng.standard_normal((N, N))
        truth[k] += PERT_STD * np.fft.ifft2(sqrt_spectrum * np.fft.fft2(noise)).real

    path = os.path.join(datadir, "sqg_truth.nc")
    write_sqg(path, ic, truth,
              "Truth for the ijedi SQG DA tests: sqg_ic.nc plus random perturbations")
    print(f"  {path}")

    # Observations of the truth on the lower boundary (level 0) and on the lid (level 1)
    for k, name in enumerate(("bottom", "lid")):
        path = os.path.join(datadir, f"obs_sqg_{name}.nc")
        write_obs(path, rng, truth[k])
        print(f"  {path}")

    # Ensemble: perturbations drawn like the truth perturbation, recentred on the background
    rng = np.random.default_rng(ENS_SEED)
    perts = np.empty((NENS, 2, N, N))
    for ie in range(NENS):
        for k in range(2):
            noise = rng.standard_normal((N, N))
            perts[ie, k] = PERT_STD * np.fft.ifft2(sqrt_spectrum * np.fft.fft2(noise)).real
    perts -= perts.mean(axis=0)
    for ie in range(NENS):
        path = os.path.join(datadir, f"ens_mem{ie+1:03d}.nc")
        write_sqg(path, ic, background + perts[ie],
                  f"Member {ie+1} of the ijedi SQG DA test ensemble")
        print(f"  {path}")
    ic.close()


if __name__ == "__main__":
    main()

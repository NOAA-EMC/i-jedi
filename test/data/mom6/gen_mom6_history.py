#!/usr/bin/env python3
"""
Generate a synthetic MOM6 history (restart) NetCDF file with realistic
ocean background fields.

Only the variables actually read by the SOCA MOM6 IO layer are written:
  - Temp   (potential temperature, degC)
  - Salt   (salinity, PPT)
  - ave_ssh (sea surface height, m)

Required dimensions (lonh, lath, Layer) and their coordinate variables are
also written so that the reader can validate grid sizes.

Usage:
    python3 gen_mom6_history.py <output_nc_path>
"""

import sys
import numpy as np
import netCDF4 as nc


# ---------------------------------------------------------------------------
# Grid definition  (must match ocean_hgrid / geometry YAML: 72x35x25)
# ---------------------------------------------------------------------------
NLON = 72   # lonh
NLAT = 35   # lath
NZ   = 25   # Layer

# h-point latitudes and longitudes (5° resolution, same as CDL)
LATH = np.array([
    -82.5, -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5,
    -32.5, -27.5, -22.5, -17.5, -12.5,  -7.5,  -2.5,   2.5,   7.5,  12.5,
     17.5,  22.5,  27.5,  32.5,  37.5,  42.5,  47.5,  52.5,  57.5,  62.5,
     67.5,  72.5,  77.5,  82.5,  87.5
])

LONH = np.array([
    -277.5, -272.5, -267.5, -262.5, -257.5, -252.5, -247.5, -242.5,
    -237.5, -232.5, -227.5, -222.5, -217.5, -212.5, -207.5, -202.5,
    -197.5, -192.5, -187.5, -182.5, -177.5, -172.5, -167.5, -162.5,
    -157.5, -152.5, -147.5, -142.5, -137.5, -132.5, -127.5, -122.5,
    -117.5, -112.5, -107.5, -102.5,  -97.5,  -92.5,  -87.5,  -82.5,
     -77.5,  -72.5,  -67.5,  -62.5,  -57.5,  -52.5,  -47.5,  -42.5,
     -37.5,  -32.5,  -27.5,  -22.5,  -17.5,  -12.5,   -7.5,   -2.5,
       2.5,    7.5,   12.5,   17.5,   22.5,   27.5,   32.5,   37.5,
      42.5,   47.5,   52.5,   57.5,   62.5,   67.5,   72.5,   77.5
])

# q-point (corner) latitudes and longitudes
LATQ = np.array([
    -80., -75., -70., -65., -60., -55., -50., -45., -40., -35.,
    -30., -25., -20., -15., -10.,  -5.,   0.,   5.,  10.,  15.,
     20.,  25.,  30.,  35.,  40.,  45.,  50.,  55.,  60.,  65.,
     70.,  75.,  80.,  85.,  90.
])

LONQ = LONH - 2.5   # staggered half-cell west of lonh

# Mid-layer depths [m] (positive down, same values as CDL)
LAYER = np.array([
      2.5,     7.5,    12.505,   17.545,   22.705,
     28.17,   34.285,  41.61,    50.99,    63.63,
     81.165, 105.755, 140.17,  187.88,   253.155,
    341.17,  458.115, 611.29,  809.23,  1061.82,
   1380.41, 1777.955, 2269.125, 2870.44, 4429.32232906817
])

# Interface depths [m] (NZ+1 values)
INTERFACE = np.array([
       0.,      5.,     10.,     15.01,   20.08,
      25.33,   31.01,   37.56,   45.66,   56.32,
      70.94,   91.39,  120.12,  160.22,  215.54,
     290.77,  391.57,  524.66,  697.92,  920.54,
    1203.1,  1557.72, 1998.19, 2540.06, 3200.82,
    5657.82465813634
])


# ---------------------------------------------------------------------------
# Synthetic ocean background fields
# ---------------------------------------------------------------------------

def make_temperature():
    """
    Potential temperature [degC].
    - Warm at the surface (~28 °C equator, ~-1.8 °C ice-covered poles)
    - Exponential thermocline with e-folding depth ~200 m
    - Uniform abyssal temperature of ~2 °C
    """
    lat_rad = np.deg2rad(LATH)

    # Surface temperature: latitudinal gradient
    T_surf = 28.0 * np.cos(lat_rad) - 1.8 * np.abs(np.sin(lat_rad))
    T_surf = np.clip(T_surf, -1.8, 30.0)   # physical limits

    T_deep = 2.0          # abyssal temperature [°C]
    z_thermo = 200.0      # thermocline e-folding depth [m]

    # Shape: (NZ, NLAT, NLON)
    T = np.empty((NZ, NLAT, NLON))
    for k, z in enumerate(LAYER):
        frac = np.exp(-z / z_thermo)
        T_level = T_surf * frac + T_deep * (1.0 - frac)   # (NLAT,)
        T[k] = T_level[:, np.newaxis]                      # broadcast over lon

    # Add a small longitudinal perturbation (mesoscale-like)
    lon_rad = np.deg2rad(LONH)
    eddy = 0.5 * np.sin(3 * lon_rad)    # ±0.5 °C wave
    T += eddy[np.newaxis, np.newaxis, :]

    return T.astype(np.float64)


def make_salinity():
    """
    Salinity [PPT / PSU].
    - Subtropical evaporation maximum (~35.8 PSU at ±20-25°)
    - Equatorial and polar freshening (~34.5 PSU)
    - Nearly uniform deep (~34.7 PSU)
    """
    lat_rad = np.deg2rad(LATH)

    # Surface salinity: double-maximum subtropical pattern
    S_surf = 35.0 + 0.8 * np.cos(2 * lat_rad) - 0.5 * np.cos(lat_rad)
    S_surf = np.clip(S_surf, 28.0, 37.0)

    S_deep = 34.7
    z_halo = 500.0   # halocline e-folding depth [m]

    S = np.empty((NZ, NLAT, NLON))
    for k, z in enumerate(LAYER):
        frac = np.exp(-z / z_halo)
        S_level = S_surf * frac + S_deep * (1.0 - frac)
        S[k] = S_level[:, np.newaxis]

    return S.astype(np.float64)


def make_ssh():
    """
    Sea surface height [m].
    - Dynamic topography: subtropical highs (~+0.5 m), equatorial/subpolar lows
    - Small-amplitude mesoscale noise (~±0.05 m)
    """
    lat_rad = np.deg2rad(LATH)
    lon_rad = np.deg2rad(LONH)

    # Background dynamic topography
    eta_bg = 0.5 * np.cos(lat_rad) ** 2 - 0.15   # (NLAT,)
    eta = eta_bg[:, np.newaxis] * np.ones((NLAT, NLON))

    # Mesoscale noise
    rng = np.random.default_rng(42)
    eta += 0.05 * rng.standard_normal((NLAT, NLON))

    return eta.astype(np.float64)


# ---------------------------------------------------------------------------
# NetCDF writer
# ---------------------------------------------------------------------------

def write(path):
    Temp   = make_temperature()   # (NZ, NLAT, NLON)
    Salt   = make_salinity()
    ave_ssh = make_ssh()           # (NLAT, NLON)

    with nc.Dataset(path, "w", format="NETCDF4_CLASSIC") as ds:
        # --- dimensions ---
        ds.createDimension("Time",      1)
        ds.createDimension("Layer",     NZ)
        ds.createDimension("Interface", NZ + 1)
        ds.createDimension("lath",      NLAT)
        ds.createDimension("lonh",      NLON)
        ds.createDimension("latq",      NLAT)
        ds.createDimension("lonq",      NLON)

        # --- coordinate variables ---
        v = ds.createVariable("Time", "f8", ("Time",))
        v.long_name = "Time"
        v.units = "days"
        v.cartesian_axis = "T"
        v[:] = [0.0]

        v = ds.createVariable("Layer", "f8", ("Layer",))
        v.long_name = "Layer z-rho"
        v.units = "meter"
        v.cartesian_axis = "Z"
        v.positive = "down"
        v[:] = LAYER

        v = ds.createVariable("Interface", "f8", ("Interface",))
        v.long_name = "Interface z-rho"
        v.units = "meter"
        v.cartesian_axis = "Z"
        v.positive = "down"
        v[:] = INTERFACE

        v = ds.createVariable("lath", "f8", ("lath",))
        v.long_name = "Latitude"
        v.units = "degrees_north"
        v.cartesian_axis = "Y"
        v[:] = LATH

        v = ds.createVariable("lonh", "f8", ("lonh",))
        v.long_name = "Longitude"
        v.units = "degrees_east"
        v.cartesian_axis = "X"
        v[:] = LONH

        v = ds.createVariable("latq", "f8", ("latq",))
        v.long_name = "Latitude"
        v.units = "degrees_north"
        v.cartesian_axis = "Y"
        v[:] = LATQ

        v = ds.createVariable("lonq", "f8", ("lonq",))
        v.long_name = "Longitude"
        v.units = "degrees_east"
        v.cartesian_axis = "X"
        v[:] = LONQ

        # --- state variables ---
        v = ds.createVariable("Temp", "f8", ("Time", "Layer", "lath", "lonh"))
        v.long_name = "Potential Temperature"
        v.units = "degC"
        v[0] = Temp

        v = ds.createVariable("Salt", "f8", ("Time", "Layer", "lath", "lonh"))
        v.long_name = "Salinity"
        v.units = "PPT"
        v[0] = Salt

        v = ds.createVariable("ave_ssh", "f8", ("Time", "lath", "lonh"))
        v.long_name = "Time average sea surface height"
        v.units = "meter"
        v[0] = ave_ssh

        # --- global attributes ---
        ds.filename = "RESTART/MOM.res.nc"
        ds.history  = "Generated by gen_mom6_history.py"

    print(f"  {path}")


# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output_nc_path>")
        sys.exit(1)
    write(sys.argv[1])


if __name__ == "__main__":
    main()

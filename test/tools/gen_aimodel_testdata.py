#!/usr/bin/env python3
"""
Generate structured latitude-longitude test data laid out the way the AI weather models
expect their inputs.

The layout follows the ARCO-ERA5 / WeatherBench-2 convention that NeuralGCM, GraphCast, Pangu
and Aurora all consume:

    dimensions   time, level, latitude, longitude
    upper air    (time, level, latitude, longitude), float32, SI units
    surface      (time, latitude, longitude), float32, SI units
    level        pressure in hPa, ascending (top of atmosphere first)
    latitude     ascending, south to north
    longitude    ascending, 0 to 360, east of Greenwich

The horizontal grid is the N32 regular Gaussian grid (128 longitudes x 64 latitudes) that the
pretrained NeuralGCM pressure-level models run on. Gaussian latitudes are computed as the
roots of the Legendre polynomial so they match what atlas builds for 'regular_gaussian N=32',
which is what makes these files readable by the ijedi 'structured netcdf' io backend.

Two files are written:

    aimodel_n32_l37.nc          the seven prognostic variables plus two surface forcings
    aimodel_n32_l37_latrev.nc   the same data with latitude descending, north to south, to
                                exercise the reader's axis-order detection

Fields are smooth analytic functions of latitude, longitude and pressure with plausible
magnitudes; they are not meteorologically meaningful, but they are reproducible, cheap and
distinguishable from one another, which is what an io test needs.

Usage: python3 gen_aimodel_testdata.py <output_dir>
"""

import os
import sys

import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    sys.exit("ERROR: netCDF4 python package is required. "
             "Install with: pip install netCDF4")

# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------
N = 32                 # Gaussian truncation: nlon = 4N = 128, nlat = 2N = 64
NLON = 4 * N
NLAT = 2 * N

# The 37 ERA5 pressure levels, hPa, top of atmosphere first.
LEVELS_HPA = np.array([
    1, 2, 3, 5, 7, 10, 20, 30, 50, 70,
    100, 125, 150, 175, 200, 225, 250, 300, 350, 400,
    450, 500, 550, 600, 650, 700, 750, 775, 800, 825,
    850, 875, 900, 925, 950, 975, 1000], dtype=np.float64)

# ERA5 names, exactly as NeuralGCM and GraphCast expect them.
UPPER_AIR = {
    "u_component_of_wind": "m s**-1",
    "v_component_of_wind": "m s**-1",
    "temperature": "K",
    "geopotential": "m**2 s**-2",
    "specific_humidity": "kg kg**-1",
    "specific_cloud_liquid_water_content": "kg kg**-1",
    "specific_cloud_ice_water_content": "kg kg**-1",
}
SURFACE = {
    "sea_surface_temperature": "K",
    "sea_ice_cover": "1",
}

VALID_TIME = "2026-01-01T00:00:00"
HOURS_SINCE_EPOCH = 491232.0   # 2026-01-01 00Z in hours since 1970-01-01


def gaussian_latitudes(nlat):
    """Gaussian latitudes in degrees, ascending (south to north).

    These are the roots of the Legendre polynomial of degree nlat, which is what a Gaussian
    quadrature grid is defined by and what atlas uses for 'regular_gaussian'.
    """
    nodes, _ = np.polynomial.legendre.leggauss(nlat)
    return np.degrees(np.arcsin(nodes))


def build_fields(lats_deg, lons_deg):
    """Smooth analytic fields on the (level, lat, lon) grid.

    Each variable gets a different combination of zonal wavenumber, meridional structure and
    vertical profile so that a transposed or level-shuffled read cannot go unnoticed.
    """
    lat = np.radians(lats_deg)[:, None]          # (nlat, 1)
    lon = np.radians(lons_deg)[None, :]          # (1, nlon)
    p = LEVELS_HPA[:, None, None] / 1000.0       # (nlev, 1, 1), sigma-like

    coslat = np.cos(lat)
    sinlat = np.sin(lat)

    fields = {}
    # Zonal jet, strongest in mid-latitudes near 250 hPa.
    fields["u_component_of_wind"] = (
        40.0 * coslat * sinlat**2 * np.exp(-((np.log(p) + 1.4) ** 2)) +
        5.0 * np.cos(2.0 * lon) * coslat) * np.ones_like(p)
    # Weaker meridional flow, shifted a quarter wavelength so u and v differ everywhere.
    fields["v_component_of_wind"] = (
        8.0 * np.sin(3.0 * lon) * coslat**2) * np.ones_like(p)
    # Temperature: warm surface tropics, cold upper-level everywhere.
    fields["temperature"] = (
        220.0 + 70.0 * p * coslat**2 + 3.0 * np.cos(lon) * coslat)
    # Geopotential from a dry hydrostatic profile, so it decreases upward like the real thing.
    fields["geopotential"] = (
        287.0 * 250.0 * np.log(1000.0 / LEVELS_HPA)[:, None, None] *
        np.ones((1, len(lats_deg), len(lons_deg))) + 500.0 * sinlat**2)
    # Moisture concentrated in the tropical lower troposphere.
    fields["specific_humidity"] = (
        0.018 * p**3 * coslat**4 * (1.0 + 0.3 * np.sin(2.0 * lon)))
    fields["specific_cloud_liquid_water_content"] = (
        2.0e-4 * p**2 * np.exp(-((lat / 0.5) ** 2)) * (1.0 + 0.5 * np.cos(4.0 * lon)))
    fields["specific_cloud_ice_water_content"] = (
        1.0e-5 * (1.0 - p) ** 2 * coslat**2 * (1.0 + 0.5 * np.sin(5.0 * lon)))

    surface = {}
    # Sea surface temperature: warm equator, freezing point at the poles.
    surface["sea_surface_temperature"] = (
        271.5 + 30.0 * coslat**3 + 1.5 * np.cos(lon) * coslat) * np.ones((len(lats_deg), 1))
    # Sea ice: a smooth ramp poleward of 60 degrees, clipped to a fraction.
    surface["sea_ice_cover"] = np.clip(
        (np.abs(lats_deg)[:, None] - 60.0) / 20.0, 0.0, 1.0) * np.ones((1, len(lons_deg)))

    return fields, surface


def write_file(path, lats_deg, fields, surface):
    """Write one file. lats_deg fixes the latitude ordering; the data is reordered to match."""
    ascending = lats_deg[0] < lats_deg[-1]
    lons_deg = np.arange(NLON) * (360.0 / NLON)

    with nc.Dataset(path, "w", format="NETCDF4") as ds:
        ds.Conventions = "CF-1.8"
        ds.title = "Synthetic AI-model input on the N32 Gaussian grid"
        ds.comment = ("Generated by gen_aimodel_testdata.py for i-jedi io tests. "
                      "Analytic fields, not a real analysis.")

        ds.createDimension("time", 1)
        ds.createDimension("level", len(LEVELS_HPA))
        ds.createDimension("latitude", NLAT)
        ds.createDimension("longitude", NLON)

        v = ds.createVariable("time", "f8", ("time",))
        v.units = "hours since 1970-01-01T00:00:00"
        v.standard_name = "time"
        v.axis = "T"
        v[:] = [HOURS_SINCE_EPOCH]

        v = ds.createVariable("level", "f8", ("level",))
        v.units = "hPa"
        v.standard_name = "air_pressure"
        v.axis = "Z"
        v[:] = LEVELS_HPA

        v = ds.createVariable("latitude", "f8", ("latitude",))
        v.units = "degrees_north"
        v.standard_name = "latitude"
        v.axis = "Y"
        v[:] = lats_deg

        v = ds.createVariable("longitude", "f8", ("longitude",))
        v.units = "degrees_east"
        v.standard_name = "longitude"
        v.axis = "X"
        v[:] = lons_deg

        # The analytic fields are built south to north; flip when writing north to south.
        def oriented(a):
            return a if ascending else a[..., ::-1, :]

        for name, units in UPPER_AIR.items():
            v = ds.createVariable(name, "f4", ("time", "level", "latitude", "longitude"),
                                  zlib=True, complevel=1)
            v.units = units
            v[0, :, :, :] = oriented(fields[name]).astype(np.float32)

        for name, units in SURFACE.items():
            v = ds.createVariable(name, "f4", ("time", "latitude", "longitude"),
                                  zlib=True, complevel=1)
            v.units = units
            v[0, :, :] = oriented(surface[name]).astype(np.float32)

    print(f"wrote {path}")


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: gen_aimodel_testdata.py <output_dir>")
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)

    lats_asc = gaussian_latitudes(NLAT)
    fields, surface = build_fields(lats_asc, np.arange(NLON) * (360.0 / NLON))

    write_file(os.path.join(outdir, "aimodel_n32_l37.nc"), lats_asc, fields, surface)
    write_file(os.path.join(outdir, "aimodel_n32_l37_latrev.nc"), lats_asc[::-1],
               fields, surface)


if __name__ == "__main__":
    main()

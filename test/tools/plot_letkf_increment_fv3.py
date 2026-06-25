#!/usr/bin/env python3

"""Plot FV3 LETKF mean increment on a global lat/lon map.

This script reads the unstructured mean background and mean analysis files
written by the I-JEDI LETKF test, computes the increment (analysis minus
background) for a chosen variable/level, interpolates it to a regular
lat/lon grid for visualization, and writes a global map.

Example
-------
python3 test/tools/plot_letkf_increment_fv3.py \
    --analysis build/ijedi/test/letkf-fv3-an.nc \
    --background build/ijedi/test/letkf-fv3-bg.nc \
    --variable air_temperature \
    --level 63 \
    --output letkf-fv3-increment-level63.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", required=True, help="Mean analysis NetCDF file")
    parser.add_argument("--background", required=True, help="Mean background NetCDF file")
    parser.add_argument("--variable", default="air_temperature",
                        help="Variable to plot from the NetCDF files")
    parser.add_argument("--level", type=int, default=0,
                        help="Zero-based vertical level to plot")
    parser.add_argument("--grid-resolution", type=float, default=1.0,
                        help="Regular lat/lon plotting grid resolution in degrees")
    parser.add_argument("--title", default=None, help="Optional custom plot title")
    parser.add_argument("--output", default=None, help="Output PNG path")
    return parser.parse_args()


def read_increment(analysis_path: Path, background_path: Path, variable: str,
                   level: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with Dataset(analysis_path) as an_ds, Dataset(background_path) as bg_ds:
        for coord in ("lon", "lat", variable):
            if coord not in an_ds.variables:
                raise KeyError(f"{analysis_path} is missing variable '{coord}'")
            if coord not in bg_ds.variables:
                raise KeyError(f"{background_path} is missing variable '{coord}'")

        lon = np.asarray(an_ds.variables["lon"][:], dtype=float)
        lat = np.asarray(an_ds.variables["lat"][:], dtype=float)
        an_var = np.asarray(an_ds.variables[variable][:], dtype=float)
        bg_var = np.asarray(bg_ds.variables[variable][:], dtype=float)

    if an_var.shape != bg_var.shape:
        raise ValueError(
            f"Analysis/background shapes differ for '{variable}': "
            f"{an_var.shape} vs {bg_var.shape}"
        )

    if an_var.ndim == 1:
        if level != 0:
            raise ValueError(f"Variable '{variable}' is 2D; only level 0 is valid")
        increment = an_var - bg_var
    elif an_var.ndim == 2:
        if level < 0 or level >= an_var.shape[1]:
            raise ValueError(
                f"Requested level {level} is outside valid range 0..{an_var.shape[1] - 1}"
            )
        increment = an_var[:, level] - bg_var[:, level]
    else:
        raise ValueError(
            f"Variable '{variable}' has unsupported rank {an_var.ndim}; expected 1 or 2"
        )

    if lon.shape != lat.shape or lon.shape[0] != increment.shape[0]:
        raise ValueError(
            "Expected lon, lat, and increment to share the same node dimension, got "
            f"{lon.shape}, {lat.shape}, {increment.shape}"
        )

    lon = np.where(lon > 180.0, lon - 360.0, lon)
    return lon, lat, increment


def interpolate_to_latlon(lon: np.ndarray, lat: np.ndarray, values: np.ndarray,
                          resolution_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon_1d = np.arange(-180.0, 180.0 + resolution_deg, resolution_deg)
    lat_1d = np.arange(-90.0, 90.0 + resolution_deg, resolution_deg)
    lon_grid, lat_grid = np.meshgrid(lon_1d, lat_1d)

    points = np.column_stack((lon, lat))
    grid = griddata(points, values, (lon_grid, lat_grid), method="linear")

    if np.isnan(grid).any():
        nearest = griddata(points, values, (lon_grid, lat_grid), method="nearest")
        grid = np.where(np.isnan(grid), nearest, grid)

    return lon_grid, lat_grid, grid


def plot_increment(lon_grid: np.ndarray, lat_grid: np.ndarray, increment_grid: np.ndarray,
                   variable: str, level: int, output_path: Path, title: str | None) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
    except ImportError as err:
        raise ImportError(
            "This plotting script requires matplotlib and cartopy in the active Python environment"
        ) from err

    absmax = np.nanpercentile(np.abs(increment_grid), 98.0)
    if not np.isfinite(absmax) or absmax == 0.0:
        absmax = 1.0

    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines(linewidth=0.6)

    mesh = ax.pcolormesh(
        lon_grid,
        lat_grid,
        increment_grid,
        transform=ccrs.PlateCarree(),
        shading="auto",
        cmap="coolwarm",
        vmin=-absmax,
        vmax=absmax,
    )

    plot_title = title or f"{variable} increment (analysis - background), level {level}"
    ax.set_title(plot_title)
    cbar = fig.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.05, shrink=0.85)
    cbar.set_label(variable)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")


def main() -> None:
    args = parse_args()
    analysis_path = Path(args.analysis)
    background_path = Path(args.background)

    if args.output is None:
        output_path = Path(f"{args.variable}_increment_level{args.level:03d}.png")
    else:
        output_path = Path(args.output)

    lon, lat, increment = read_increment(
        analysis_path, background_path, args.variable, args.level
    )
    lon_grid, lat_grid, increment_grid = interpolate_to_latlon(
        lon, lat, increment, args.grid_resolution
    )
    plot_increment(
        lon_grid, lat_grid, increment_grid, args.variable, args.level, output_path, args.title
    )


if __name__ == "__main__":
    main()

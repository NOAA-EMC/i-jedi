#!/usr/bin/env python3
"""Stage NeuralGCM initial conditions and a reference forecast for the i-jedi tests.

Writes two netCDF files in the structured latitude-longitude layout that i-jedi's
'structured netcdf' io backend reads (dimensions time, level, latitude, longitude; CF
coordinate variables; SI units; ERA5 variable names):

    neuralgcm_ic_tl63.nc    initial condition on NeuralGCM's own TL63 grid
    neuralgcm_ref_tl63.nc   forecast from that initial condition, produced entirely in
                            Python by NeuralGCM itself

The reference file is the point of the exercise. Running a forecast through i-jedi and
finding that it completes proves only that the plumbing does not crash; comparing it against
a rollout that never went near C++ is what shows the bridge preserves the model's answer.

By default the initial condition comes from the ERA5 snapshot packaged inside the neuralgcm
wheel, conservatively regridded to TL63, so nothing is downloaded and the result is
reproducible. That snapshot is coarse (TL31 upsampled) and is not a good forecast initial
condition scientifically - it is a good *test* initial condition, which is a different thing.
Pass --era5 to read a real analysis instead once one is staged.

The grid is not a coincidence: NeuralGCM's TL63 grid and atlas's 'regular_gaussian N=32'
agree to 3e-14 degrees, which is why i-jedi can run this model on its own native grid with
no interpolation anywhere.

Usage: stage_neuralgcm_data.py <output_dir> [--hours 6] [--seed 0] [--era5 FILE]
"""

import argparse
import os
import sys

import numpy as np

try:
    import jax
    import neuralgcm
    import xarray
    from neuralgcm import demo
except ImportError as exc:  # pragma: no cover - environment problem, not a code path
    sys.exit(f"ERROR: the NeuralGCM environment is not importable ({exc}).\n"
             "Create it with tools/setup_python_model_env.sh and run this with that "
             "venv's python.")

# The seven prognostic variables and two surface forcings the pretrained pressure-level
# models exchange, under ERA5's names.
UPPER_AIR = [
    "u_component_of_wind",
    "v_component_of_wind",
    "temperature",
    "geopotential",
    "specific_humidity",
    "specific_cloud_liquid_water_content",
    "specific_cloud_ice_water_content",
]
SURFACE = ["sea_surface_temperature", "sea_ice_cover"]


def to_ijedi_layout(ds):
    """Reorder to (time, level, latitude, longitude) with CF coordinate attributes.

    NeuralGCM works in (time, level, longitude, latitude) - longitude-major, the transpose
    of the ERA5 convention - so this is a real reordering, not a relabelling.
    """
    out = xarray.Dataset()
    for name in list(ds.data_vars):
        var = ds[name]
        # NeuralGCM's output carries bookkeeping variables such as sim_time that have no
        # horizontal dimensions. They are not fields and are not written.
        if "latitude" not in var.dims or "longitude" not in var.dims:
            continue
        dims = ("time", "level", "latitude", "longitude") if "level" in var.dims \
            else ("time", "latitude", "longitude")
        out[name] = var.transpose(*dims)
    out = out.assign_coords({k: v for k, v in ds.coords.items() if k in
                             ("time", "level", "latitude", "longitude")})

    out.latitude.attrs.update(units="degrees_north", standard_name="latitude", axis="Y")
    out.longitude.attrs.update(units="degrees_east", standard_name="longitude", axis="X")
    if "level" in out.coords:
        out.level.attrs.update(units="hPa", standard_name="air_pressure", axis="Z")
    out.attrs["Conventions"] = "CF-1.8"
    out.attrs["comment"] = ("Staged for i-jedi tests by stage_neuralgcm_data.py; "
                            "NeuralGCM TL63 grid (= atlas regular_gaussian N=32)")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir")
    parser.add_argument("--hours", type=int, default=6,
                        help="length of the reference forecast (default 6)")
    parser.add_argument("--seed", type=int, default=0,
                        help="rng seed for the stochastic model; must match the yaml")
    parser.add_argument("--era5", default=None,
                        help="read the initial condition from this file instead of the "
                             "packaged demo snapshot")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("loading the packaged TL63 checkpoint...")
    model = neuralgcm.PressureLevelModel.from_checkpoint(
        demo.load_checkpoint_tl63_stochastic())
    coords = model.data_coords
    print(f"  internal timestep {model.timestep}, "
          f"{coords.vertical.layers} levels, "
          f"{len(coords.horizontal.longitudes)}x{len(coords.horizontal.latitudes)} grid")

    if args.era5:
        print(f"reading initial condition from {args.era5}")
        ds = xarray.load_dataset(args.era5)
        if "time" not in ds.dims:
            ds = ds.expand_dims("time")
    else:
        print("regridding the packaged ERA5 snapshot to TL63...")
        ds = demo.load_data(coords)

    keep = [v for v in UPPER_AIR + SURFACE if v in ds.data_vars]
    missing = set(UPPER_AIR + SURFACE) - set(keep)
    if missing:
        sys.exit(f"ERROR: the source data is missing {sorted(missing)}")
    ds = ds[keep]

    ic_path = os.path.join(args.outdir, "neuralgcm_ic_tl63.nc")
    to_ijedi_layout(ds).to_netcdf(ic_path)
    print(f"wrote {ic_path}")

    # --- reference forecast, entirely in Python -------------------------------------------
    print(f"running a {args.hours} h reference forecast (seed {args.seed})...")
    # encode takes forcings at a single time; unroll takes a time series, because it
    # interpolates the forcings to each internal step. Passing the single-time object to
    # unroll gives 0-d arrays and a confusing failure inside jnp.interp.
    inputs = model.inputs_from_xarray(ds.isel(time=0))
    forcings_now = model.forcings_from_xarray(ds.isel(time=0))
    forcings_series = model.forcings_from_xarray(ds)
    rng = jax.random.key(args.seed)
    state = model.encode(inputs, forcings_now, rng)

    steps = int(args.hours * 3600 // int(model.timestep.total_seconds())) \
        if hasattr(model.timestep, "total_seconds") else \
        int(np.timedelta64(args.hours, "h") / model.timestep)
    print(f"  {steps} internal steps")

    # Deliberately the same call the adapter makes - one output at the end of the whole
    # interval, rather than one per internal step. The two are mathematically equivalent but
    # give XLA different scan structures and so different floating-point summation orders,
    # which shows up as a ~1e-5 relative difference in float32. Matching the call isolates
    # the bridge as the only thing that differs between this reference and the i-jedi run,
    # which is what the reference is for.
    final_state, predictions = model.unroll(
        state, forcings_series, steps=1,
        timedelta=np.timedelta64(int(args.hours) * 3600, "s"),
        start_with_input=False)

    times = np.array([ds.time.values[0] + np.timedelta64(int(args.hours) * 3600, "s")])
    ref = model.data_to_xarray(predictions, times=times).isel(time=-1).expand_dims("time")
    # Carry the forcings through so the reference file holds the same variables as the
    # forecast i-jedi writes.
    for name in SURFACE:
        if name in ds:
            ref[name] = ds[name].isel(time=0).expand_dims("time")

    ref_path = os.path.join(args.outdir, "neuralgcm_ref_tl63.nc")
    to_ijedi_layout(ref).to_netcdf(ref_path)
    print(f"wrote {ref_path}")

    for name in UPPER_AIR[:3]:
        a, b = np.asarray(ds[name].isel(time=0)), np.asarray(ref[name].isel(time=0))
        print(f"  {name:38s} initial rms {np.sqrt((a**2).mean()):12.5e} "
              f"-> forecast rms {np.sqrt((b**2).mean()):12.5e}")


if __name__ == "__main__":
    main()

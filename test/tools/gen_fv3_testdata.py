#!/usr/bin/env python3
"""
Generate FV3 test data.

Outputs:
    1) FV3 history files: atmf006.nc and sfcf006.nc
    2) FV3/FMS restart files for a single tile domain:
         coupler.res, fv_core.res.nc, fv_core.res.tile1.nc,
         fv_tracer.res.tile1.nc, fv_srf_wnd.res.tile1.nc

The random seed is fixed so output is reproducible across runs.

Usage: python3 gen_fv3_testdata.py <output_dir>
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
# Grid / vertical coordinate parameters
# ---------------------------------------------------------------------------
NPX = 13          # grid_xt = npx - 1 = 12
NPY = 13          # grid_yt = npy - 1 = 12
NPZ = 127         # number of full levels
NTILES = 6
NX = NPX - 1      # 12
NY = NPY - 1      # 12
NMEMS = 3
HALO = 3

# GFS L127 ak/bk (interface pressures: p_half = ak + bk * ps)
AK = np.array([
    9.990000e-01, 1.605000e+00, 2.532000e+00, 3.924000e+00,
    5.976000e+00, 8.947000e+00, 1.317700e+01, 1.909600e+01,
    2.724300e+01, 3.827600e+01, 5.298400e+01, 7.229300e+01,
    9.726900e+01, 1.291100e+02, 1.691350e+02, 2.187670e+02,
    2.795060e+02, 3.528940e+02, 4.404810e+02, 5.437820e+02,
    6.642360e+02, 8.031640e+02, 9.617340e+02, 1.140931e+03,
    1.341538e+03, 1.564119e+03, 1.809028e+03, 2.076415e+03,
    2.366252e+03, 2.678372e+03, 3.012510e+03, 3.368363e+03,
    3.745646e+03, 4.144164e+03, 4.563881e+03, 5.004995e+03,
    5.468017e+03, 5.953848e+03, 6.463864e+03, 7.000000e+03,
    7.563494e+03, 8.150661e+03, 8.756529e+03, 9.376141e+03,
    1.000455e+04, 1.063685e+04, 1.126816e+04, 1.189364e+04,
    1.250852e+04, 1.310809e+04, 1.368773e+04, 1.424289e+04,
    1.476915e+04, 1.526220e+04, 1.571786e+04, 1.613209e+04,
    1.650102e+04, 1.682094e+04, 1.708832e+04, 1.729985e+04,
    1.745308e+04, 1.754835e+04, 1.758677e+04, 1.756970e+04,
    1.749870e+04, 1.737556e+04, 1.720230e+04, 1.698114e+04,
    1.671450e+04, 1.640502e+04, 1.605549e+04, 1.566886e+04,
    1.524825e+04, 1.479687e+04, 1.431804e+04, 1.381515e+04,
    1.329163e+04, 1.275092e+04, 1.219647e+04, 1.163166e+04,
    1.105983e+04, 1.048421e+04, 9.907927e+03, 9.333967e+03,
    8.765155e+03, 8.204142e+03, 7.653387e+03, 7.115147e+03,
    6.591468e+03, 6.084176e+03, 5.594876e+03, 5.124949e+03,
    4.675554e+03, 4.247633e+03, 3.841918e+03, 3.458933e+03,
    3.099010e+03, 2.762297e+03, 2.448768e+03, 2.158238e+03,
    1.890375e+03, 1.644712e+03, 1.420661e+03, 1.217528e+03,
    1.034524e+03, 8.707780e+02, 7.253480e+02, 5.972350e+02,
    4.853920e+02, 3.887340e+02, 3.061490e+02, 2.365020e+02,
    1.786510e+02, 1.314470e+02, 9.374000e+01, 6.439200e+01,
    4.227400e+01, 2.627400e+01, 1.530200e+01, 8.287000e+00,
    4.190000e+00, 1.994000e+00, 8.100000e-01, 2.320000e-01,
    2.900000e-02, 0.000000e+00, 0.000000e+00, 0.000000e+00,
], dtype=np.float32)

BK = np.array([
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    1.018000e-05, 8.141000e-05, 2.746900e-04, 6.507800e-04,
    1.270090e-03, 2.192480e-03, 3.477130e-03, 5.182280e-03,
    7.365040e-03, 1.008120e-02, 1.338492e-02, 1.732857e-02,
    2.196239e-02, 2.733428e-02, 3.348954e-02, 4.047056e-02,
    4.831661e-02, 5.706358e-02, 6.674372e-02, 7.738548e-02,
    8.900629e-02, 1.015940e-01, 1.151262e-01, 1.295762e-01,
    1.449129e-01, 1.611008e-01, 1.780999e-01, 1.958660e-01,
    2.143511e-01, 2.335031e-01, 2.532663e-01, 2.735822e-01,
    2.943890e-01, 3.156229e-01, 3.372180e-01, 3.591072e-01,
    3.812224e-01, 4.034951e-01, 4.258572e-01, 4.482413e-01,
    4.705813e-01, 4.928130e-01, 5.148743e-01, 5.367062e-01,
    5.582525e-01, 5.794605e-01, 6.002815e-01, 6.206707e-01,
    6.405875e-01, 6.599957e-01, 6.788633e-01, 6.971631e-01,
    7.148720e-01, 7.319713e-01, 7.484465e-01, 7.642871e-01,
    7.794867e-01, 7.940422e-01, 8.079541e-01, 8.212263e-01,
    8.338652e-01, 8.458801e-01, 8.572826e-01, 8.680866e-01,
    8.783077e-01, 8.879632e-01, 8.970718e-01, 9.056532e-01,
    9.137284e-01, 9.213187e-01, 9.284464e-01, 9.351338e-01,
    9.414037e-01, 9.472789e-01, 9.527821e-01, 9.579360e-01,
    9.627630e-01, 9.672851e-01, 9.715240e-01, 9.755009e-01,
    9.792364e-01, 9.827508e-01, 9.860625e-01, 9.891851e-01,
    9.921299e-01, 9.949077e-01, 9.975282e-01, 1.000000e+00,
], dtype=np.float32)

# Reference pressure levels (mb) for pfull and phalf coordinate variables
PHALF_REF = np.linspace(0.0, 1000.0, NPZ + 1)
PFULL_REF = 0.5 * (PHALF_REF[:-1] + PHALF_REF[1:])


def _axis_values(size):
    """Create FMS-like axis coordinates (1, 11, 21, ...)."""
    return (1.0 + 10.0 * np.arange(size, dtype=np.float32))


def _set_time_and_axis_attrs(vtime):
    vtime.cartesian_axis = "T"
    vtime.units = "time level"
    vtime.long_name = "Time"


def _write_restart_coupler(filepath):
    """Write a coupler.res file matching FV3 restart text layout."""
    with open(filepath, "w", encoding="ascii") as fobj:
        fobj.write("     2        (Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)\n")
        fobj.write("  2024     5    26    23     0     0        Model start time:   year, month, day, hour, minute, second\n")
        fobj.write("  2024     5    27     0     0     0        Current model time: year, month, day, hour, minute, second\n")


def generate_fv_core_res(filepath):
    """Generate fv_core.res.nc containing vertical coordinate arrays."""
    ds = nc.Dataset(filepath, "w", format="NETCDF4")

    nlevp1 = NPZ + 1
    ds.createDimension("Time", 1)
    ds.createDimension("xaxis_1", nlevp1)

    v = ds.createVariable("xaxis_1", "f4", ("xaxis_1",), fill_value=np.nan)
    v.axis = "X"
    v[:] = _axis_values(nlevp1)

    v = ds.createVariable("Time", "f4", ("Time",), fill_value=np.nan)
    _set_time_and_axis_attrs(v)
    v[:] = np.array([1.0], dtype=np.float32)

    v = ds.createVariable("ak", "f4", ("Time", "xaxis_1"), fill_value=np.nan)
    v[:] = np.asarray(AK[:nlevp1], dtype=np.float32)[np.newaxis, :]

    v = ds.createVariable("bk", "f4", ("Time", "xaxis_1"), fill_value=np.nan)
    v[:] = np.asarray(BK[:nlevp1], dtype=np.float32)[np.newaxis, :]

    ds.close()


def generate_fv_core_res_tile1(filepath, rng):
    """Generate fv_core.res.tile1.nc for a one-tile restart domain."""
    ds = nc.Dataset(filepath, "w", format="NETCDF4")

    # Restart files contain only the compute domain (no halos).
    # FMS domain IO reads/writes compute data and checks dimension sizes
    # against the global compute extent (sum across all tiles).
    x1 = NX        # npx-1 = 12  cell-centred x
    y2 = NY        # npy-1 = 12  cell-centred y
    x2 = NX + 1   # npx   = 13  D-grid v (extra x)
    y1 = NY + 1   # npy   = 13  D-grid u (extra y)

    ds.createDimension("Time", None)
    ds.createDimension("zaxis_1", NPZ)
    ds.createDimension("yaxis_2", y2)
    ds.createDimension("xaxis_1", x1)
    ds.createDimension("yaxis_1", y1)
    ds.createDimension("xaxis_2", x2)

    ds.setncattr("NumFilesInSet", np.int32(1))

    vx1 = ds.createVariable("xaxis_1", "f4", ("xaxis_1",), fill_value=np.nan)
    vx1.axis = "X"
    vx1[:] = _axis_values(x1)

    vx2 = ds.createVariable("xaxis_2", "f4", ("xaxis_2",), fill_value=np.nan)
    vx2.axis = "X"
    vx2[:] = _axis_values(x2)

    vy1 = ds.createVariable("yaxis_1", "f4", ("yaxis_1",), fill_value=np.nan)
    vy1.axis = "Y"
    vy1[:] = _axis_values(y1)

    vy2 = ds.createVariable("yaxis_2", "f4", ("yaxis_2",), fill_value=np.nan)
    vy2.axis = "Y"
    vy2[:] = _axis_values(y2)

    vz = ds.createVariable("zaxis_1", "f4", ("zaxis_1",), fill_value=np.nan)
    vz.axis = "Z"
    vz[:] = np.arange(1, NPZ + 1, dtype=np.float32)

    vt = ds.createVariable("Time", "f4", ("Time",), fill_value=np.nan)
    _set_time_and_axis_attrs(vt)
    vt[:] = np.array([1.0], dtype=np.float32)

    def add3d(name, dims, mean, std, vmin=None, vmax=None):
        var = ds.createVariable(name, "f4", dims, fill_value=np.nan)
        shp = tuple(len(ds.dimensions[d]) for d in dims)
        data = rng.normal(mean, std, shp).astype(np.float32)
        if vmin is not None:
            data = np.clip(data, vmin, None)
        if vmax is not None:
            data = np.clip(data, None, vmax)
        var[:] = data

    add3d("DZ", ("Time", "zaxis_1", "yaxis_2", "xaxis_1"), -615.0, 500.0, vmin=-3500.0, vmax=-10.0)
    add3d("T", ("Time", "zaxis_1", "yaxis_2", "xaxis_1"), 248.0, 30.0, vmin=160.0, vmax=330.0)
    add3d("W", ("Time", "zaxis_1", "yaxis_2", "xaxis_1"), 0.0, 1.5)
    add3d("DELP", ("Time", "zaxis_1", "yaxis_2", "xaxis_1"), 778.0, 300.0, vmin=0.1)
    add3d("phis", ("Time", "yaxis_2", "xaxis_1"), 219.0, 600.0, vmin=-50.0, vmax=6000.0)
    add3d("u", ("Time", "zaxis_1", "yaxis_1", "xaxis_1"), 7.0, 15.0)
    add3d("ua", ("Time", "zaxis_1", "yaxis_2", "xaxis_1"), 7.0, 15.0)
    add3d("v", ("Time", "zaxis_1", "yaxis_2", "xaxis_2"), 0.0, 15.0)
    add3d("va", ("Time", "zaxis_1", "yaxis_2", "xaxis_1"), 0.0, 15.0)

    ds.close()


def generate_fv_tracer_res_tile1(filepath, rng):
    """Generate fv_tracer.res.tile1.nc with standard tracer fields."""
    ds = nc.Dataset(filepath, "w", format="NETCDF4")

    # Compute domain only – no halo.
    x1 = NX        # 12
    y1 = NY        # 12

    ds.createDimension("Time", None)
    ds.createDimension("xaxis_1", x1)
    ds.createDimension("yaxis_1", y1)
    ds.createDimension("zaxis_1", NPZ)

    ds.setncattr("NumFilesInSet", np.int32(1))

    vx1 = ds.createVariable("xaxis_1", "f4", ("xaxis_1",), fill_value=np.nan)
    vx1.axis = "X"
    vx1[:] = _axis_values(x1)

    vy1 = ds.createVariable("yaxis_1", "f4", ("yaxis_1",), fill_value=np.nan)
    vy1.axis = "Y"
    vy1[:] = _axis_values(y1)

    vz = ds.createVariable("zaxis_1", "f4", ("zaxis_1",), fill_value=np.nan)
    vz.axis = "Z"
    vz[:] = np.arange(1, NPZ + 1, dtype=np.float32)

    vt = ds.createVariable("Time", "f4", ("Time",), fill_value=np.nan)
    _set_time_and_axis_attrs(vt)
    vt[:] = np.array([1.0], dtype=np.float32)

    tracer_fields = [
        "sphum", "liq_wat", "ice_wat", "rainwat", "snowwat", "graupel",
        "water_nc", "ice_nc", "rain_nc", "o3mr", "liq_aero", "ice_aero",
        "sgs_tke", "smoke", "dust", "coarsepm",
    ]

    for name in tracer_fields:
        var = ds.createVariable(name, "f4", ("Time", "zaxis_1", "yaxis_1", "xaxis_1"), fill_value=np.nan)
        if name == "o3mr":
            data = rng.normal(2.1e-6, 3.0e-6, (1, NPZ, y1, x1)).astype(np.float32)
            data = np.clip(data, 1e-8, 2e-5)
        elif name == "sphum":
            data = rng.normal(0.0024, 0.003, (1, NPZ, y1, x1)).astype(np.float32)
            data = np.clip(data, 1e-8, 0.025)
        else:
            data = rng.normal(0.0, 1.0e-7, (1, NPZ, y1, x1)).astype(np.float32)
            data = np.clip(data, 0.0, None)
        var[:] = data

    ds.close()


def generate_fv_srf_wnd_res_tile1(filepath, rng):
    """Generate fv_srf_wnd.res.tile1.nc with surface wind fields."""
    ds = nc.Dataset(filepath, "w", format="NETCDF4")

    # Compute domain only – no halo.
    x1 = NX        # 12
    y1 = NY        # 12

    ds.createDimension("Time", None)
    ds.createDimension("xaxis_1", x1)
    ds.createDimension("yaxis_1", y1)

    ds.setncattr("NumFilesInSet", np.int32(1))

    vx1 = ds.createVariable("xaxis_1", "f4", ("xaxis_1",), fill_value=np.nan)
    vx1.axis = "X"
    vx1[:] = _axis_values(x1)

    vy1 = ds.createVariable("yaxis_1", "f4", ("yaxis_1",), fill_value=np.nan)
    vy1.axis = "Y"
    vy1[:] = _axis_values(y1)

    vt = ds.createVariable("Time", "f4", ("Time",), fill_value=np.nan)
    _set_time_and_axis_attrs(vt)
    vt[:] = np.array([1.0], dtype=np.float32)

    us = ds.createVariable("u_srf", "f4", ("Time", "yaxis_1", "xaxis_1"), fill_value=np.nan)
    us[:] = rng.normal(7.0, 8.0, (1, y1, x1)).astype(np.float32)

    vs = ds.createVariable("v_srf", "f4", ("Time", "yaxis_1", "xaxis_1"), fill_value=np.nan)
    vs[:] = rng.normal(0.0, 8.0, (1, y1, x1)).astype(np.float32)

    ds.close()


def _write_fixed_string_char_array(var, text):
    """Write a Python string into a fixed-width netCDF char variable."""
    n = var.shape[-1]
    # Keep trailing bytes as NUL padding rather than explicit spaces.
    var[...] = nc.stringtochar(np.array([text[:n]], dtype=f"S{n}"))


def _axis_centers(n):
    return np.arange(1, n + 1, dtype=np.float64)


def _axis_interfaces(n):
    return np.arange(1, n + 1, dtype=np.float64)


def generate_input_tile_grid(filepath, rng):
    """Generate INPUT/C775_grid.tile7.halo3.nc style grid file."""
    # FMS grid readers are most robust with classic NetCDF encodings.
    ds = nc.Dataset(filepath, "w", format="NETCDF3_64BIT_OFFSET")

    nx_core = NX
    ny_core = NY
    # FV3 read_grid reads supergrid corner coordinates for the compute grid
    # plus halo. For a regional tile this is 2 * (n + 2*halo) + 1 points.
    nxp = 2 * (nx_core + 2 * HALO) + 1
    nyp = 2 * (ny_core + 2 * HALO) + 1
    nx_h = nxp - 1
    ny_h = nyp - 1

    ds.createDimension("string", 255)
    ds.createDimension("nx", nx_h)
    ds.createDimension("ny", ny_h)
    ds.createDimension("nxp", nxp)
    ds.createDimension("nyp", nyp)

    ds.setncattr("history", "gen_fv3_testdata.py gnomonic_ed")
    ds.setncattr("source", "FV3GFS")
    ds.setncattr("grid", "akappa")
    ds.setncattr("plat", np.float64(38.5))
    ds.setncattr("plon", np.float64(-97.5))
    ds.setncattr("pazi", np.float64(0.0))
    ds.setncattr("delx", np.float64(0.05796706665290621))
    ds.setncattr("dely", np.float64(0.05827354522162335))
    ds.setncattr("lx", np.int32(-(nx_h // 2)))
    ds.setncattr("ly", np.int32(-(ny_h // 2)))
    ds.setncattr("a", np.float64(0.1208865850407068))
    ds.setncattr("k", np.float64(-0.3416958316117138))
    ds.setncattr("avg_cell_size", np.float64(12906.134028879402))
    ds.setncattr("RES_equiv", np.int32(775))

    tile = ds.createVariable("tile", "S1", ("string",))
    tile.standard_name = "grid_tile_spec"
    tile[:] = nc.stringtochar(np.array(["tile7".ljust(255)], dtype="S255"))[0, :]

    lon0 = 236.0
    lat0 = 20.0
    dlon = 0.058
    dlat = 0.058
    rot_deg = 13.0

    xif = _axis_interfaces(nxp)
    yif = _axis_interfaces(nyp)
    xgrid, ygrid = np.meshgrid(xif, yif)
    xloc = xgrid - np.mean(xif)
    yloc = ygrid - np.mean(yif)

    rot = np.deg2rad(rot_deg)
    xr = np.cos(rot) * xloc - np.sin(rot) * yloc
    yr = np.sin(rot) * xloc + np.cos(rot) * yloc

    lon = lon0 + dlon * xr + 2.5e-4 * xloc * yloc
    lat = lat0 + dlat * yr + 1.8e-4 * (xloc ** 2 - yloc ** 2)
    lon += rng.normal(0.0, 7.0e-4, lon.shape)
    lat += rng.normal(0.0, 7.0e-4, lat.shape)
    lon = np.mod(lon, 360.0)
    lat = np.clip(lat, -89.9, 89.9)

    x = ds.createVariable("x", "f8", ("nyp", "nxp"))
    x.standard_name = "geographic_longitude"
    x.units = "degree_east"
    x.hstagger = "C"
    x[:] = lon

    y = ds.createVariable("y", "f8", ("nyp", "nxp"))
    y.standard_name = "geographic_latitude"
    y.units = "degree_north"
    y.hstagger = "C"
    y[:] = lat

    # Compute consistent metrics from the generated lon/lat corner grid.
    earth_radius = 6_371_000.0

    lon_rad = np.deg2rad(lon)
    lat_rad = np.deg2rad(lat)

    dlon_x = np.diff(lon_rad, axis=1)
    dlat_x = np.diff(lat_rad, axis=1)
    dlon_x = (dlon_x + np.pi) % (2.0 * np.pi) - np.pi
    lat_mid_x = 0.5 * (lat_rad[:, 1:] + lat_rad[:, :-1])
    dx_vals = earth_radius * np.hypot(dlat_x, np.cos(lat_mid_x) * dlon_x)

    dlon_y = np.diff(lon_rad, axis=0)
    dlat_y = np.diff(lat_rad, axis=0)
    dlon_y = (dlon_y + np.pi) % (2.0 * np.pi) - np.pi
    lat_mid_y = 0.5 * (lat_rad[1:, :] + lat_rad[:-1, :])
    dy_vals = earth_radius * np.hypot(dlat_y, np.cos(lat_mid_y) * dlon_y)

    dx_center = 0.5 * (dx_vals[:-1, :] + dx_vals[1:, :])
    dy_center = 0.5 * (dy_vals[:, :-1] + dy_vals[:, 1:])
    area_vals = dx_center * dy_center

    area = ds.createVariable("area", "f8", ("ny", "nx"))
    area.standard_name = "grid_cell_area"
    area.units = "m2_area"
    area.hstagger = "H"
    area[:] = area_vals

    dx = ds.createVariable("dx", "f8", ("nyp", "nx"))
    dx.standard_name = "dx"
    dx.units = "meters"
    dx.hstagger = "H"
    dx[:] = dx_vals

    dy = ds.createVariable("dy", "f8", ("ny", "nxp"))
    dy.standard_name = "dy"
    dy.units = "meters"
    dy.hstagger = "H"
    dy[:] = dy_vals

    # Orientation of local axes relative to east in degrees.
    dlonx_full = np.empty_like(lon_rad)
    dlatx_full = np.empty_like(lat_rad)
    dlonx_full[:, :-1] = lon_rad[:, 1:] - lon_rad[:, :-1]
    dlatx_full[:, :-1] = lat_rad[:, 1:] - lat_rad[:, :-1]
    dlonx_full[:, -1] = lon_rad[:, -1] - lon_rad[:, -2]
    dlatx_full[:, -1] = lat_rad[:, -1] - lat_rad[:, -2]
    dlonx_full = (dlonx_full + np.pi) % (2.0 * np.pi) - np.pi
    lat_for_xang = lat_rad
    angle_dx_vals = np.rad2deg(np.arctan2(dlatx_full, np.cos(lat_for_xang) * dlonx_full))

    dlony_full = np.empty_like(lon_rad)
    dlaty_full = np.empty_like(lat_rad)
    dlony_full[:-1, :] = lon_rad[1:, :] - lon_rad[:-1, :]
    dlaty_full[:-1, :] = lat_rad[1:, :] - lat_rad[:-1, :]
    dlony_full[-1, :] = lon_rad[-1, :] - lon_rad[-2, :]
    dlaty_full[-1, :] = lat_rad[-1, :] - lat_rad[-2, :]
    dlony_full = (dlony_full + np.pi) % (2.0 * np.pi) - np.pi
    angle_dy_vals = np.rad2deg(np.arctan2(dlaty_full, np.cos(lat_rad) * dlony_full)) - 90.0

    angle_dx = ds.createVariable("angle_dx", "f8", ("nyp", "nxp"))
    angle_dx.standard_name = "angle_dx"
    angle_dx.units = "degrees"
    angle_dx.hstagger = "C"
    angle_dx[:] = angle_dx_vals

    angle_dy = ds.createVariable("angle_dy", "f8", ("nyp", "nxp"))
    angle_dy.standard_name = "angle_dy"
    angle_dy.units = "degrees"
    angle_dy.hstagger = "C"
    angle_dy[:] = angle_dy_vals

    ds.setncattr("delx", np.float64(np.mean(np.abs(np.diff(lon, axis=1)))))
    ds.setncattr("dely", np.float64(np.mean(np.abs(np.diff(lat, axis=0)))))
    ds.setncattr("avg_cell_size", np.float64(np.sqrt(np.mean(area_vals))))

    ds.close()


def generate_input_grid_spec(filepath):
    """Generate INPUT/grid_spec.nc style metadata for one tile."""
    ds = nc.Dataset(filepath, "w", format="NETCDF3_64BIT_OFFSET")

    ds.createDimension("ntiles", 1)
    ds.createDimension("string", 255)

    ds.setncattr("grid_version", "0.2")
    ds.setncattr("code_version", "$Name: fre-nctools-bronx-10 $")
    ds.setncattr("history", "gen_fv3_testdata.py")

    mosaic = ds.createVariable("mosaic", "S1", ("string",))
    mosaic.standard_name = "grid_mosaic_spec"
    mosaic.children = "gridtiles"
    mosaic.contact_regions = "contacts"
    mosaic.grid_descriptor = ""

    gridlocation = ds.createVariable("gridlocation", "S1", ("string",))
    gridlocation.standard_name = "grid_file_location"

    gridfiles = ds.createVariable("gridfiles", "S1", ("ntiles", "string"))
    gridtiles = ds.createVariable("gridtiles", "S1", ("ntiles", "string"))

    _write_fixed_string_char_array(mosaic, "C775_mosaic.halo3")
    _write_fixed_string_char_array(gridlocation, "./")

    n = len(ds.dimensions["string"])
    gridfiles[0, :] = nc.stringtochar(np.array(["C775_grid.tile7.halo3.nc"[:n]], dtype=f"S{n}"))[0, :]
    gridtiles[0, :] = nc.stringtochar(np.array(["tile7"[:n]], dtype=f"S{n}"))[0, :]

    ds.close()


def _write_global_attrs(ds):
    """Write global attributes common to both atm and sfc files."""
    ds.setncattr("ak", AK)
    ds.setncattr("bk", BK)
    ds.source = "FV3-JEDI"
    ds.grid = "cubed_sphere"
    ds.grid_id = np.int32(1)
    ds.ncnsto = np.int32(9)
    ds.hydrostatic = "non-hydrostatic"


def _write_common_dims_and_coords(ds, rng):
    """Create dimensions and coordinate variables shared by atm/sfc files."""
    ds.createDimension("grid_xt", NX)
    ds.createDimension("grid_yt", NY)
    ds.createDimension("pfull", NPZ)
    ds.createDimension("phalf", NPZ + 1)
    ds.createDimension("tile", NTILES)
    ds.createDimension("time", 1)
    ds.createDimension("nchars", 20)

    # tile
    v = ds.createVariable("tile", "i4", ("tile",))
    v.long_name = "cubed-sphere face"
    v[:] = np.arange(1, NTILES + 1, dtype=np.int32)

    # grid_xt
    v = ds.createVariable("grid_xt", "f8", ("grid_xt",))
    v.cartesian_axis = "X"
    v[:] = np.arange(1, NX + 1, dtype=np.float64)

    # grid_yt
    v = ds.createVariable("grid_yt", "f8", ("grid_yt",))
    v.cartesian_axis = "Y"
    v[:] = np.arange(1, NY + 1, dtype=np.float64)

    # lon (degrees_E, 0..360)
    v = ds.createVariable("lon", "f8", ("tile", "grid_yt", "grid_xt"))
    v.long_name = "T-cell longitude"
    v.units = "degrees_E"
    lon_base = np.linspace(1.0, 358.0, NX)
    lon_2d = np.broadcast_to(lon_base[np.newaxis, :], (NY, NX)).copy()
    for t in range(NTILES):
        offset = t * 60.0  # rotate each tile
        v[t, :, :] = (lon_2d + offset + rng.uniform(-0.5, 0.5, (NY, NX))) % 360.0

    # lat (degrees_N, -90..90)
    v = ds.createVariable("lat", "f8", ("tile", "grid_yt", "grid_xt"))
    v.long_name = "T-cell latitude"
    v.units = "degrees_N"
    lat_base = np.linspace(-84.0, 84.0, NY)
    lat_2d = np.broadcast_to(lat_base[:, np.newaxis], (NY, NX)).copy()
    for t in range(NTILES):
        v[t, :, :] = lat_2d + rng.uniform(-0.5, 0.5, (NY, NX))

    # pfull
    v = ds.createVariable("pfull", "f8", ("pfull",))
    v.positive = "down"
    v.long_name = "ref full pressure level"
    v.units = "mb"
    v.edges = "phalf"
    v.cartesian_axis = "Z"
    v[:] = PFULL_REF

    # phalf
    v = ds.createVariable("phalf", "f8", ("phalf",))
    v.positive = "down"
    v.long_name = "ref half pressure level"
    v.units = "mb"
    v.cartesian_axis = "Z"
    v[:] = PHALF_REF

    # time
    v = ds.createVariable("time", "f8", ("time",))
    v.long_name = "time"
    v.calendar = "JULIAN"
    v.calendar_type = "JULIAN"
    v.cartesian_axis = "T"
    v.units = "hours since 2021-03-22 12:00:00"
    v[:] = [0.0]

    # time_iso
    v = ds.createVariable("time_iso", "S1", ("time", "nchars"))
    v.long_name = "valid time"
    v.description = "ISO 8601 datetime string"
    iso_str = "2020-12-15T00:00:00Z"
    v[0, :] = np.array(list(iso_str), dtype="S1")


def _add_3d_field(ds, name, long_name, units, standard_name, rng,
                  mean, std, vmin=None, vmax=None):
    """Add a 3-D atmospheric field (time, tile, pfull, grid_yt, grid_xt)."""
    v = ds.createVariable(name, "f8",
                          ("time", "tile", "pfull", "grid_yt", "grid_xt"))
    v.long_name = long_name
    v.units = units
    v.standard_name = standard_name
    v.coordinates = "lon lat"
    v.grid_mapping = "cubed_sphere"
    data = rng.normal(mean, std, (1, NTILES, NPZ, NY, NX))
    if vmin is not None:
        data = np.clip(data, vmin, None)
    if vmax is not None:
        data = np.clip(data, None, vmax)
    v[:] = data


def _add_2d_field(ds, name, long_name, units, standard_name, rng,
                  mean, std, vmin=None, vmax=None):
    """Add a 2-D atmospheric field (time, tile, grid_yt, grid_xt)."""
    v = ds.createVariable(name, "f8",
                          ("time", "tile", "grid_yt", "grid_xt"))
    v.long_name = long_name
    v.units = units
    v.standard_name = standard_name
    v.coordinates = "lon lat"
    v.grid_mapping = "cubed_sphere"
    data = rng.normal(mean, std, (1, NTILES, NY, NX))
    if vmin is not None:
        data = np.clip(data, vmin, None)
    if vmax is not None:
        data = np.clip(data, None, vmax)
    v[:] = data


def generate_atmf(filepath, rng):
    """Generate the atmosphere history file."""
    ds = nc.Dataset(filepath, "w", format="NETCDF4")

    _write_global_attrs(ds)
    _write_common_dims_and_coords(ds, rng)

    # ---- 3-D fields (time, tile, pfull, grid_yt, grid_xt) ----

    # eastward_wind (m/s): roughly -60..80, mean ~7
    _add_3d_field(ds, "ugrd", "eastward_wind", "ms-1", "eastward_wind",
                  rng, mean=7.0, std=15.0)

    # northward_wind (m/s): roughly -75..75, mean ~0
    _add_3d_field(ds, "vgrd", "northward_wind", "ms-1", "northward_wind",
                  rng, mean=0.0, std=15.0)

    # air_temperature (K): 173..313, mean ~248
    _add_3d_field(ds, "tmp", "air_temperature", "K", "air_temperature",
                  rng, mean=248.0, std=30.0, vmin=160.0, vmax=330.0)

    # specific_humidity (kg/kg): ~0..0.02, mean ~0.0024
    _add_3d_field(ds, "spfh", "specific_humidity", "kgkg-1",
                  "specific_humidity",
                  rng, mean=0.0024, std=0.003, vmin=1e-7, vmax=0.025)

    # layer_thickness (m): negative values, ~-3000..-16, mean ~-615
    _add_3d_field(ds, "delz", "layer_thickness", "m", "layer_thickness",
                  rng, mean=-615.0, std=500.0, vmin=-3500.0, vmax=-10.0)

    # air_pressure_thickness (Pa): ~0.6..1771, mean ~778
    _add_3d_field(ds, "dpres", "air_pressure_thickness", "pa",
                  "air_pressure_thickness",
                  rng, mean=778.0, std=300.0, vmin=0.1)

    # ---- 2-D fields (time, tile, grid_yt, grid_xt) ----

    # surface_pressure (Pa): ~54000..104000, mean ~98776
    _add_2d_field(ds, "pressfc", "surface_pressure", "Pa",
                  "surface_pressure",
                  rng, mean=98776.0, std=5000.0, vmin=40000.0, vmax=110000.0)

    # sfc_geopotential_height (m): ~-27..5100, mean ~219
    _add_2d_field(ds, "hgtsfc",
                  "sfc_geopotential_height_times_grav", "m",
                  "sfc_geopotential_height_times_grav",
                  rng, mean=219.0, std=600.0, vmin=-50.0, vmax=6000.0)

    # ozone_mass_mixing_ratio (kg/kg): ~3e-8..1.8e-5, mean ~2.1e-6
    _add_3d_field(ds, "o3mr", "ozone_mass_mixing_ratio", "kgkg-1",
                  "ozone_mass_mixing_ratio",
                  rng, mean=2.1e-6, std=3.0e-6, vmin=1e-8, vmax=2e-5)

    ds.close()


def generate_sfcf(filepath, rng):
    """Generate the surface history file (dimensions + global attrs only)."""
    ds = nc.Dataset(filepath, "w", format="NETCDF4")
    _write_global_attrs(ds)

    # Same dimensions as the atm file, but no data variables
    ds.createDimension("grid_xt", NX)
    ds.createDimension("grid_yt", NY)
    ds.createDimension("pfull", NPZ)
    ds.createDimension("phalf", NPZ + 1)
    ds.createDimension("tile", NTILES)
    ds.createDimension("time", 1)
    ds.createDimension("nchars", 20)

    ds.close()


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output_dir>")
        sys.exit(1)

    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)

    # Generate test data for deterministic (imem == 0) and ensemble members (imem > 0)
    seed = 20201215
    for imem in range(0, NMEMS+1):
        rng = np.random.default_rng(seed=seed)
        suffix = ''
        if imem > 0:
            memstr = str(imem).zfill(3)
            suffix = f'_mem{memstr}'

        atm_path = os.path.join(outdir, f"atmf006{suffix}.nc")
        generate_atmf(atm_path, rng)
        print(f"  {atm_path}")

        sfc_path = os.path.join(outdir, f"sfcf006{suffix}.nc")
        generate_sfcf(sfc_path, rng)
        print(f"  {sfc_path}")

        core_path = os.path.join(outdir, f"fv_core.res{suffix}.nc")
        generate_fv_core_res(core_path)
        print(f"  {core_path}")

        core_tile_path = os.path.join(outdir, f"fv_core.res.tile1{suffix}.nc")
        generate_fv_core_res_tile1(core_tile_path, rng)
        print(f"  {core_tile_path}")

        trcr_tile_path = os.path.join(outdir, f"fv_tracer.res.tile1{suffix}.nc")
        generate_fv_tracer_res_tile1(trcr_tile_path, rng)
        print(f"  {trcr_tile_path}")

        sfcw_tile_path = os.path.join(outdir, f"fv_srf_wnd.res.tile1{suffix}.nc")
        generate_fv_srf_wnd_res_tile1(sfcw_tile_path, rng)
        print(f"  {sfcw_tile_path}")
        
        seed += 1

    coupler_path = os.path.join(outdir, f"coupler.res")
    _write_restart_coupler(coupler_path)
    print(f"  {coupler_path}")

    # Write FV3 grid metadata into a top-level INPUT directory
    # (sibling of data_generated/) so FMS can find it via INPUT paths.
    input_dir = os.path.abspath(os.path.join(outdir, "..", "..", "INPUT"))
    os.makedirs(input_dir, exist_ok=True)

    input_tile_path = os.path.join(input_dir, "C775_grid.tile7.halo3.nc")
    generate_input_tile_grid(input_tile_path, rng)
    print(f"  {input_tile_path}")

    input_spec_path = os.path.join(input_dir, "grid_spec.nc")
    generate_input_grid_spec(input_spec_path)
    print(f"  {input_spec_path}")

    print(f"FV3 test data written to: {outdir}")


if __name__ == "__main__":
    main()

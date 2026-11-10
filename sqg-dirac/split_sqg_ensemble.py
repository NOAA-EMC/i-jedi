#!/usr/bin/env python3
"""
Split an sqgturb EnKF ensemble file into per-member SQG model files readable by i-jedi IoSQG.

Usage: split_sqg_ensemble.py <ensemble file> <output directory> [variable]

The ensemble file has the background (or analysis) boundary PV of each member,
<variable>(t, ens, z, y, x) in K (default variable: pv_b), at a single time (the last record is
used). IoSQG reads pv(t, z, y, x) in m/s and converts it to K with the variable's scale_to_K
attribute, so each member is written as
    <output directory>/ens/sqg_ens_memNNN.nc  (NNN = 001, 002, ...)
with pv = <variable> / (f*theta0/g), as written by the SQG model; the ensemble mean is written
the same way to <output directory>/sqg_ens_mean.nc.
"""

import os
import sys

import netCDF4 as nc
import numpy as np


def write_sqg(path, src, pv, scale, comment):
    ds = nc.Dataset(path, "w", format="NETCDF4")
    ds.createDimension("t", None)
    for name, size in zip(("z", "y", "x"), pv.shape):
        ds.createDimension(name, size)
    v = ds.createVariable("pv", "f4", ("t", "z", "y", "x"))
    v.units = "m/s"
    v.scale_to_K = np.float32(scale)
    v[0] = (pv / scale).astype(np.float32)
    ds.setncatts({a: src.getncattr(a) for a in src.ncattrs() if a not in ("history", "NCO")})
    ds.comment = comment
    ds.close()


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    infile, outdir = sys.argv[1], sys.argv[2]
    varname = sys.argv[3] if len(sys.argv) == 4 else "pv_b"

    src = nc.Dataset(infile)
    scale = float(src.f) * float(src.theta0) / float(src.g)
    ens = np.asarray(src[varname][-1], dtype=np.float64)  # (ens, z, y, x) in K
    nens = ens.shape[0]

    os.makedirs(os.path.join(outdir, "ens"), exist_ok=True)
    for ie in range(nens):
        path = os.path.join(outdir, "ens", f"sqg_ens_mem{ie+1:03d}.nc")
        write_sqg(path, src, ens[ie], scale,
                  f"member {ie+1} of {varname} in {os.path.basename(infile)}")
    path = os.path.join(outdir, "sqg_ens_mean.nc")
    write_sqg(path, src, ens.mean(axis=0), scale,
              f"ensemble mean of {varname} ({nens} members) in {os.path.basename(infile)}")
    src.close()
    print(f"Wrote {nens} members and the ensemble mean of {varname} to {outdir} "
          f"(scale_to_K = {scale:.6g})")


if __name__ == "__main__":
    main()

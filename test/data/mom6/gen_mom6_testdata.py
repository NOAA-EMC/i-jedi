#!/usr/bin/env python3
"""
Generate MOM6 geometry test data from CDL files.

Recreates ocean_hgrid.nc and ocean_topog.nc from CDL dumps of the real
72x35x25 GFDL low-resolution tripolar grid.

CDL files are committed to the repo (no git-lfs required). ncgen recreates
the netCDF files at test time.  Grid parameters (NIGLOBAL, NJGLOBAL, etc.)
are specified directly in the geometry YAML config — no MOM_input file is
generated or needed.

Usage: python3 gen_mom6_testdata.py <output_dir>
"""

import os
import sys
import subprocess

def ncgen(cdl_path, nc_path):
    result = subprocess.run(
        ["ncgen", "-o", nc_path, cdl_path],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ncgen failed for {cdl_path}:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output_dir>")
        sys.exit(1)

    outdir = sys.argv[1]
    input_dir = os.path.join(outdir, "INPUT")
    os.makedirs(input_dir, exist_ok=True)

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # -------------------------------------------------------------------------
    # ocean_hgrid.nc and ocean_topog.nc from CDL
    # -------------------------------------------------------------------------
    for name in ("ocean_hgrid", "ocean_topog"):
        cdl = os.path.join(script_dir, f"{name}.cdl")
        nc  = os.path.join(input_dir, f"{name}.nc")
        ncgen(cdl, nc)
        print(f"  {nc}")

    # -------------------------------------------------------------------------
    # SOCA reference gridspec (used by comparison test)
    # -------------------------------------------------------------------------
    soca_cdl = os.path.join(script_dir, "soca_gridspec.72x35x25.cdl")
    soca_nc  = os.path.join(outdir, "soca_gridspec.72x35x25.nc")
    ncgen(soca_cdl, soca_nc)
    print(f"  {soca_nc}")

    print(f"MOM6 test data written to: {outdir}")


if __name__ == "__main__":
    main()

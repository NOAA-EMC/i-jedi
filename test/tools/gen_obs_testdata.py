#!/usr/bin/env python3
"""
Generate observation test data from CDL files.

CDL files are committed to the repo (no git-lfs required). ncgen recreates
the netCDF files at test time.

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
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_dir> <output_dir>")
        sys.exit(1)

    indir = sys.argv[1]
    outdir = sys.argv[2]
    print("Generating observation test data in:", outdir)
    os.makedirs(outdir, exist_ok=True)

    # -------------------------------------------------------------------------
    # observations
    # -------------------------------------------------------------------------
    for name in ("sst", "aircraft", "insitu_temp_profile_argo_3prof"):
        cdl = os.path.join(indir, f"{name}.cdl")
        nc  = os.path.join(outdir, f"{name}.nc")
        ncgen(cdl, nc)
        print(f"  {nc}")

if __name__ == "__main__":
    main()

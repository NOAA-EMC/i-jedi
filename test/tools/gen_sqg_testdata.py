#!/usr/bin/env python3
"""
Generate SQG model test data from CDL files.

sqg_ic.cdl is a 64x64 SQG state (sqgturb-cpp layout: pv(t,z,y,x) in m/s with a
scale_to_K attribute), committed as CDL so no binary files are needed; ncgen
recreates the netCDF file at test time. The CDL keeps 9 significant digits, so
the float32 values round-trip exactly.

Usage: python3 gen_sqg_testdata.py <input_dir> <output_dir>
"""

import os
import subprocess
import sys


def main():
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    indir, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    for name in ["sqg_ic"]:
        result = subprocess.run(
            ["ncgen", "-o", os.path.join(outdir, name + ".nc"),
             os.path.join(indir, name + ".cdl")],
            capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ncgen failed for {name}.cdl:\n{result.stderr}", file=sys.stderr)
            sys.exit(1)
    print(f"SQG test data written to {outdir}")


if __name__ == "__main__":
    main()

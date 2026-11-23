#!/usr/bin/env python3
"""
Generate WW3 test data: a small unstructured mesh in the WW3 gmsh layout and
2-D analytic fields on its nodes.

Nodes: a 4-degree lon-lat grid (lon in 0-360, 82S-86N, offset from the MOM6
test grid), triangulated two triangles per cell, periodic in longitude (so
triangles cross the 0/360 seam). Nodes over land of the 72x35 MOM6 test
topography are dropped with their triangles, which gives coastlines.

Writes, in <mom6 data dir>:
  ww3_test.msh   gmsh v2 ASCII, as fix/wave/<grid>.msh: "id lon lat depth",
                 elements "id 2 3 0 1 0 n1 n2 n3"
  ww3_fields.nc  hs, ice on dimension nx (msh node order)

Usage: gen_ww3_testdata.py <data_generated/mom6> <data_generated/ww3>
"""

import os
import sys

import netCDF4
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_mom6_testdata import ice_concentration  # noqa: E402


def wave_height(lon, lat):
    """Smooth analytic significant wave height (m)."""
    return 2.0 + np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon))


def xyz(lon, lat):
    lam, phi = np.deg2rad(lon), np.deg2rad(lat)
    return np.stack([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)], -1)


def main():
    mom6_dir, out_dir = sys.argv[1], sys.argv[2]

    # MOM6 test grid T-points and depth (land = 0)
    with netCDF4.Dataset(os.path.join(mom6_dir, "soca_gridspec.72x35x25.nc")) as gs:
        mlon = gs["lon"][0].data.ravel()
        mlat = gs["lat"][0].data.ravel()
    with netCDF4.Dataset(os.path.join(mom6_dir, "ocean_topog.nc")) as ds:
        mdepth = ds["depth"][:].data.ravel()

    # Lon-lat node grid
    lons = np.arange(2.0, 360.0, 4.0)
    lats = np.arange(-82.0, 87.0, 4.0)
    ni, nj = lons.size, lats.size
    LON, LAT = np.meshgrid(lons, lats)
    lon, lat = LON.ravel(), LAT.ravel()

    # Depth of the nearest MOM6 T-point; nodes over land are dropped
    p, q = xyz(lon, lat), xyz(mlon, mlat)
    depth = mdepth[np.argmax(p @ q.T, axis=1)]
    wet = depth > 0.0

    # Triangles of each (periodic) cell, kept when all nodes are wet
    j, i = np.meshgrid(np.arange(nj - 1), np.arange(ni), indexing="ij")
    sw, se = j * ni + i, j * ni + (i + 1) % ni
    nw, ne = sw + ni, se + ni
    tri = np.concatenate([np.stack([sw, se, ne], -1).reshape(-1, 3),
                          np.stack([sw, ne, nw], -1).reshape(-1, 3)])
    tri = tri[wet[tri].all(axis=1)]

    # Keep the nodes of the kept triangles, renumbered 1..N
    keep = np.zeros(lon.size, bool)
    keep[tri] = True
    newid = np.cumsum(keep)
    tri = newid[tri]
    lon, lat, depth = lon[keep], lat[keep], depth[keep]

    with open(os.path.join(out_dir, "ww3_test.msh"), "w") as f:
        f.write("$MeshFormat\n2 0 8\n$EndMeshFormat\n")
        f.write(f"$Nodes\n{lon.size}\n")
        for n in range(lon.size):
            f.write(f"{n + 1} {lon[n]:.6f} {lat[n]:.6f} {depth[n]:.6f}\n")
        f.write(f"$EndNodes\n$Elements\n{len(tri)}\n")
        for e, (a, b, c) in enumerate(tri):
            f.write(f"{e + 1} 2 3 0 1 0 {a} {b} {c}\n")
        f.write("$EndElements\n")

    with netCDF4.Dataset(os.path.join(out_dir, "ww3_fields.nc"), "w") as ds:
        ds.createDimension("time", 1)
        ds.createDimension("ny", 1)
        ds.createDimension("nx", lon.size)
        ds.createVariable("hs", "f8", ("time", "nx"))[:] = wave_height(lon, lat)[None]
        ds.createVariable("ice", "f4", ("time", "ny", "nx"),
                          fill_value=-999.0)[:] = ice_concentration(lat)[None, None]

    print(f"ww3_test.msh: {lon.size} nodes, {len(tri)} triangles")


if __name__ == "__main__":
    main()

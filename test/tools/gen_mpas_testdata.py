#!/usr/bin/env python3
"""
Generate MPAS test data: small MPAS Voronoi meshes, built from a subdivided icosahedron,
with the mesh variables of an MPAS invariant/static file and graph partition files.

Outputs:
    1) x1.642.invariant.nc          global quasi-uniform mesh (642 cells)
       x1.642.graph.info.part.N     partitions for N = 2, 6
    2) regional.invariant.nc        limited-area cut of the global mesh, with bdyMask fields
       regional.graph.info.part.N   partitions for N = 2, 6

The meshes are deterministic so output is reproducible across runs.

Usage: python3 gen_mpas_testdata.py <output_dir>
"""

import os
import sys
import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    sys.exit("ERROR: netCDF4 python package is required. "
             "Install with: pip install netCDF4")

SPHERE_RADIUS = 6371229.0   # MPAS sphere radius (m)
N_SUBDIVISIONS = 3          # 10 * 4**3 + 2 = 642 cells
MAX_EDGES = 10
N_VERT_LEVELS = 6
Z_TOP = 30000.0
PARTITIONS = (2, 6)
# Limited-area cut: cells within this angular radius of the centre
REGIONAL_CENTER_LATLON = (40.0, -100.0)
REGIONAL_RADIUS_DEG = 55.0
N_BDY_LAYERS = 7


# ---------------------------------------------------------------------------
# Global icosahedral Voronoi mesh
# ---------------------------------------------------------------------------
def icosahedron():
    t = (1.0 + 5.0 ** 0.5) / 2.0
    pts = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
                    [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
                    [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], dtype=float)
    faces = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
             [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
             [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
             [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]]
    return pts / np.linalg.norm(pts, axis=1)[:, None], faces


def subdivide(pts, faces):
    pts = list(pts)
    midpoint = {}

    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in midpoint:
            p = pts[a] + pts[b]
            pts.append(p / np.linalg.norm(p))
            midpoint[key] = len(pts) - 1
        return midpoint[key]

    new_faces = []
    for a, b, c in faces:
        ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
        new_faces += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
    return np.array(pts), new_faces


def spherical_triangle_area(a, b, c):
    # Unit-sphere triangle area (L'Huilier via the atan2 formula)
    num = np.abs(np.einsum("ij,ij->i", a, np.cross(b, c)))
    den = 1.0 + np.einsum("ij,ij->i", a, b) + np.einsum("ij,ij->i", b, c) \
        + np.einsum("ij,ij->i", c, a)
    return 2.0 * np.arctan2(num, den)


def global_mesh():
    pts, faces = icosahedron()
    for _ in range(N_SUBDIVISIONS):
        pts, faces = subdivide(pts, faces)
    faces = np.array(faces)

    # Orient every triangle counter-clockwise seen from outside the sphere
    a, b, c = pts[faces[:, 0]], pts[faces[:, 1]], pts[faces[:, 2]]
    flip = np.einsum("ij,ij->i", np.cross(b - a, c - a), a) < 0
    faces[flip] = faces[flip][:, [0, 2, 1]]

    # Voronoi vertices = triangle circumcentres
    a, b, c = pts[faces[:, 0]], pts[faces[:, 1]], pts[faces[:, 2]]
    vert = np.cross(b - a, c - a)
    vert /= np.linalg.norm(vert, axis=1)[:, None]

    n_cells = len(pts)
    tris_of_cell = [[] for _ in range(n_cells)]
    for v, f in enumerate(faces):
        for cell in f:
            tris_of_cell[cell].append(v)

    vertices_on_cell = np.zeros((n_cells, MAX_EDGES), dtype=np.int32)
    cells_on_cell = np.zeros((n_cells, MAX_EDGES), dtype=np.int32)
    n_edges_on_cell = np.zeros(n_cells, dtype=np.int32)
    area_cell = np.zeros(n_cells)
    for cell, tris in enumerate(tris_of_cell):
        # Order the surrounding Voronoi vertices counter-clockwise around the cell centre
        p = pts[cell]
        e1 = np.cross(p, [0.0, 0.0, 1.0] if abs(p[2]) < 0.9 else [1.0, 0.0, 0.0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(p, e1)
        ang = [np.arctan2(vert[t] @ e2, vert[t] @ e1) for t in tris]
        tris = [tris[i] for i in np.argsort(ang)]
        n = len(tris)
        n_edges_on_cell[cell] = n
        vertices_on_cell[cell, :n] = np.array(tris) + 1
        # Cell across the edge between consecutive vertices
        for j in range(n):
            shared = set(faces[tris[j]]) & set(faces[tris[(j + 1) % n]])
            shared.discard(cell)
            cells_on_cell[cell, j] = shared.pop() + 1
        nxt = [tris[(j + 1) % n] for j in range(n)]
        area_cell[cell] = spherical_triangle_area(
            np.repeat(p[None, :], n, axis=0), vert[tris], vert[nxt]).sum()

    lat_cell = np.arcsin(np.clip(pts[:, 2], -1.0, 1.0))
    lon_cell = np.mod(np.arctan2(pts[:, 1], pts[:, 0]), 2.0 * np.pi)
    lat_vertex = np.arcsin(np.clip(vert[:, 2], -1.0, 1.0))
    lon_vertex = np.mod(np.arctan2(vert[:, 1], vert[:, 0]), 2.0 * np.pi)

    return dict(latCell=lat_cell, lonCell=lon_cell, areaCell=area_cell * SPHERE_RADIUS ** 2,
                latVertex=lat_vertex, lonVertex=lon_vertex,
                nEdgesOnCell=n_edges_on_cell, cellsOnCell=cells_on_cell,
                verticesOnCell=vertices_on_cell, cellsOnVertex=(faces + 1).astype(np.int32))


# ---------------------------------------------------------------------------
# Limited-area cut, renumbered as the MPAS-Limited-Area tool does (0 = missing)
# ---------------------------------------------------------------------------
def regional_mesh(mesh):
    lat0, lon0 = np.radians(REGIONAL_CENTER_LATLON)
    cosd = (np.sin(lat0) * np.sin(mesh["latCell"]) +
            np.cos(lat0) * np.cos(mesh["latCell"]) * np.cos(mesh["lonCell"] - lon0))
    keep_cell = np.degrees(np.arccos(np.clip(cosd, -1, 1))) < REGIONAL_RADIUS_DEG
    keep_vertex = np.zeros(len(mesh["latVertex"]), dtype=bool)
    keep_vertex[mesh["verticesOnCell"][keep_cell][mesh["verticesOnCell"][keep_cell] > 0] - 1] \
        = True

    cell_map = np.zeros(len(keep_cell) + 1, dtype=np.int32)
    cell_map[1:][keep_cell] = np.arange(1, keep_cell.sum() + 1)
    vertex_map = np.zeros(len(keep_vertex) + 1, dtype=np.int32)
    vertex_map[1:][keep_vertex] = np.arange(1, keep_vertex.sum() + 1)

    reg = {k: mesh[k][keep_cell] for k in ("latCell", "lonCell", "areaCell", "nEdgesOnCell")}
    reg.update({k: mesh[k][keep_vertex] for k in ("latVertex", "lonVertex")})
    reg["cellsOnCell"] = cell_map[mesh["cellsOnCell"][keep_cell]]
    reg["verticesOnCell"] = vertex_map[mesh["verticesOnCell"][keep_cell]]
    reg["cellsOnVertex"] = cell_map[mesh["cellsOnVertex"][keep_vertex]]

    # bdyMaskCell: 7 at the outermost cells, decreasing inwards, 0 in the interior
    n_cells = keep_cell.sum()
    bdy_cell = np.zeros(n_cells, dtype=np.int32)
    layer = np.array([(row[:n] == 0).any() for row, n in
                      zip(reg["cellsOnCell"], reg["nEdgesOnCell"])])
    for value in range(N_BDY_LAYERS, 0, -1):
        bdy_cell[layer & (bdy_cell == 0)] = value
        nbrs = reg["cellsOnCell"][layer]
        nxt = np.zeros(n_cells, dtype=bool)
        nxt[nbrs[nbrs > 0] - 1] = True
        layer = nxt & (bdy_cell == 0)
    reg["bdyMaskCell"] = bdy_cell
    # bdyMaskVertex: 7 where a cell of the vertex lies outside the domain
    reg["bdyMaskVertex"] = np.where((reg["cellsOnVertex"] == 0).any(axis=1), 7,
                                    0).astype(np.int32)
    return reg


# ---------------------------------------------------------------------------
# Vertical grid and output
# ---------------------------------------------------------------------------
def zgrid(lat, lon):
    # Smooth terrain and terrain-following levels decaying to a flat model top
    terrain = 1500.0 * np.maximum(0.0, np.cos(2.0 * lat) * np.sin(3.0 * lon))
    eta = np.linspace(0.0, 1.0, N_VERT_LEVELS + 1)
    return terrain[:, None] * (1.0 - eta[None, :]) + Z_TOP * eta[None, :]


def partition(lon, nparts):
    # Contiguous longitude wedges of (nearly) equal size
    owner = np.empty(len(lon), dtype=int)
    owner[np.argsort(lon, kind="stable")] = np.arange(len(lon)) * nparts // len(lon)
    return owner


def write_mesh(path, mesh):
    n_cells, n_vertices = len(mesh["latCell"]), len(mesh["latVertex"])
    with nc.Dataset(path, "w", format="NETCDF4") as f:
        f.createDimension("Time", None)
        f.createDimension("nCells", n_cells)
        f.createDimension("nVertices", n_vertices)
        f.createDimension("maxEdges", MAX_EDGES)
        f.createDimension("vertexDegree", 3)
        f.createDimension("nVertLevels", N_VERT_LEVELS)
        f.createDimension("nVertLevelsP1", N_VERT_LEVELS + 1)
        f.on_sphere = "YES"
        f.sphere_radius = SPHERE_RADIUS
        f.is_periodic = "NO"

        def var(name, dtype, dims, data, units):
            v = f.createVariable(name, dtype, dims)
            v.units = units
            v[:] = data

        var("latCell", "f4", ("nCells",), mesh["latCell"], "rad")
        var("lonCell", "f4", ("nCells",), mesh["lonCell"], "rad")
        var("areaCell", "f4", ("nCells",), mesh["areaCell"], "m^2")
        var("indexToCellID", "i4", ("nCells",), np.arange(1, n_cells + 1), "-")
        var("nEdgesOnCell", "i4", ("nCells",), mesh["nEdgesOnCell"], "-")
        var("cellsOnCell", "i4", ("nCells", "maxEdges"), mesh["cellsOnCell"], "-")
        var("verticesOnCell", "i4", ("nCells", "maxEdges"), mesh["verticesOnCell"], "-")
        var("latVertex", "f4", ("nVertices",), mesh["latVertex"], "rad")
        var("lonVertex", "f4", ("nVertices",), mesh["lonVertex"], "rad")
        var("indexToVertexID", "i4", ("nVertices",), np.arange(1, n_vertices + 1), "-")
        var("cellsOnVertex", "i4", ("nVertices", "vertexDegree"), mesh["cellsOnVertex"], "-")
        var("zgrid", "f4", ("nCells", "nVertLevelsP1"),
            zgrid(mesh["latCell"], mesh["lonCell"]), "m")
        for name in ("bdyMaskCell", "bdyMaskVertex"):
            if name in mesh:
                dim = "nCells" if name.endswith("Cell") else "nVertices"
                var(name, "i4", (dim,), mesh[name], "-")


def write_partitions(prefix, lon):
    for nparts in PARTITIONS:
        np.savetxt(f"{prefix}{nparts}", partition(lon, nparts), fmt="%d")


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output_dir>")
        sys.exit(1)
    outdir = sys.argv[1]
    print("Generating MPAS test data in:", outdir)
    os.makedirs(outdir, exist_ok=True)

    mesh = global_mesh()
    n = len(mesh["latCell"])
    write_mesh(os.path.join(outdir, f"x1.{n}.invariant.nc"), mesh)
    write_partitions(os.path.join(outdir, f"x1.{n}.graph.info.part."), mesh["lonCell"])

    reg = regional_mesh(mesh)
    write_mesh(os.path.join(outdir, "regional.invariant.nc"), reg)
    write_partitions(os.path.join(outdir, "regional.graph.info.part."), reg["lonCell"])
    print(f"  global: {n} cells, regional: {len(reg['latCell'])} cells")


if __name__ == "__main__":
    main()

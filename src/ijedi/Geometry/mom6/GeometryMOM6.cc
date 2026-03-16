// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include <iomanip>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <netcdf.h>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Metadata.h"

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

#include "ijedi/Geometry/mom6/GeometryMOM6.h"

// -------------------------------------------------------------------------------------------------
namespace ijedi
{

// -------------------------------------------------------------------------------------------------
// Replicates FMS compute_extent(): integer-division partition with remainder
// distributed to lower-rank processes.
int GeometryMOM6::computeExtent(int N, int ndivs, int pe)
{
  int base  = N / ndivs;
  int extra = N % ndivs;
  return base + (pe < extra ? 1 : 0);
}

// -------------------------------------------------------------------------------------------------
// Read T-grid lon/lat from ocean_hgrid.nc — full global array, size [njGlobal * niGlobal].
// The supergrid is (2*NJ+1) × (2*NI+1); T-centres are at stride-2 odd indices.
void GeometryMOM6::readHgridLonLat(const std::string & path,
                                    std::vector<double> & lon,
                                    std::vector<double> & lat) const
{
  int ncid;
  if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
    throw eckit::CantOpenFile(path, Here());

  const int nxp = 2 * niGlobal_ + 1;
  const int nyp = 2 * njGlobal_ + 1;

  std::vector<double> xFull(static_cast<size_t>(nyp) * nxp);
  std::vector<double> yFull(static_cast<size_t>(nyp) * nxp);

  int xid, yid;
  nc_inq_varid(ncid, "x", &xid);
  nc_inq_varid(ncid, "y", &yid);
  nc_get_var_double(ncid, xid, xFull.data());
  nc_get_var_double(ncid, yid, yFull.data());
  nc_close(ncid);

  lon.resize(static_cast<size_t>(njGlobal_) * niGlobal_);
  lat.resize(static_cast<size_t>(njGlobal_) * niGlobal_);

  for (int jG = 0; jG < njGlobal_; ++jG) {
    const int jSuper = 2 * jG + 1;
    for (int iG = 0; iG < niGlobal_; ++iG) {
      const int iSuper = 2 * iG + 1;
      lon[jG * niGlobal_ + iG] = xFull[jSuper * nxp + iSuper];
      lat[jG * niGlobal_ + iG] = yFull[jSuper * nxp + iSuper];
    }
  }
}

// -------------------------------------------------------------------------------------------------
// Read T-cell metrics from ocean_hgrid.nc — full global arrays, size [njGlobal * niGlobal].
//
// Supergrid variable dimensions and T-cell formulas (0-based iG, jG):
//   dx  (nyp, nx=2*NI)   : dxT  = dx[2*jG+1, 2*iG] + dx[2*jG+1, 2*iG+1]
//   dy  (ny=2*NJ, nxp)   : dyT  = dy[2*jG,   2*iG+1] + dy[2*jG+1, 2*iG+1]
//   area(ny=2*NJ, nx=2*NI): areaT = sum of 2×2 block area[2*jG:2*jG+2, 2*iG:2*iG+2]
void GeometryMOM6::readHgridMetrics(const std::string & path,
                                     std::vector<double> & dxT,
                                     std::vector<double> & dyT,
                                     std::vector<double> & areaT) const
{
  int ncid;
  if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
    throw eckit::CantOpenFile(path, Here());

  const int nx   = 2 * niGlobal_;
  const int nyp  = 2 * njGlobal_ + 1;
  const int nxp  = 2 * niGlobal_ + 1;
  const int nySub = 2 * njGlobal_;

  std::vector<double> dxFull  (static_cast<size_t>(nyp)  * nx);
  std::vector<double> dyFull  (static_cast<size_t>(nySub) * nxp);
  std::vector<double> areaFull(static_cast<size_t>(nySub) * nx);

  int vid;
  nc_inq_varid(ncid, "dx",   &vid); nc_get_var_double(ncid, vid, dxFull.data());
  nc_inq_varid(ncid, "dy",   &vid); nc_get_var_double(ncid, vid, dyFull.data());
  nc_inq_varid(ncid, "area", &vid); nc_get_var_double(ncid, vid, areaFull.data());
  nc_close(ncid);

  dxT.resize(static_cast<size_t>(njGlobal_) * niGlobal_);
  dyT.resize(static_cast<size_t>(njGlobal_) * niGlobal_);
  areaT.resize(static_cast<size_t>(njGlobal_) * niGlobal_);

  for (int jG = 0; jG < njGlobal_; ++jG) {
    const int j2 = 2 * jG;
    for (int iG = 0; iG < niGlobal_; ++iG) {
      const int i2 = 2 * iG;
      const int n  = jG * niGlobal_ + iG;

      dxT[n]   = dxFull[(j2 + 1) * nx + i2] + dxFull[(j2 + 1) * nx + i2 + 1];
      dyT[n]   = dyFull[j2 * nxp + (i2 + 1)] + dyFull[(j2 + 1) * nxp + (i2 + 1)];
      areaT[n] = areaFull[ j2      * nx + i2]
               + areaFull[ j2      * nx + i2 + 1]
               + areaFull[(j2 + 1) * nx + i2]
               + areaFull[(j2 + 1) * nx + i2 + 1];
    }
  }
}

// -------------------------------------------------------------------------------------------------
// Read U/V-cell centre lon/lat from ocean_hgrid.nc — full global arrays, size [njGlobal * niGlobal].
//
// MOM6 C-grid staggering (0-based iG, jG):
//   U-point (east face of T-cell): x/y[2*jG+1, 2*iG+2]
//   V-point (north face of T-cell): x/y[2*jG+2, 2*iG+1]
void GeometryMOM6::readHgridUVLonLat(const std::string & path,
                                      std::vector<double> & lonU,
                                      std::vector<double> & latU,
                                      std::vector<double> & lonV,
                                      std::vector<double> & latV) const
{
  int ncid;
  if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
    throw eckit::CantOpenFile(path, Here());

  const int nxp = 2 * niGlobal_ + 1;
  const int nyp = 2 * njGlobal_ + 1;

  std::vector<double> xFull(static_cast<size_t>(nyp) * nxp);
  std::vector<double> yFull(static_cast<size_t>(nyp) * nxp);

  int xid, yid;
  nc_inq_varid(ncid, "x", &xid);
  nc_inq_varid(ncid, "y", &yid);
  nc_get_var_double(ncid, xid, xFull.data());
  nc_get_var_double(ncid, yid, yFull.data());
  nc_close(ncid);

  lonU.resize(static_cast<size_t>(njGlobal_) * niGlobal_);
  latU.resize(static_cast<size_t>(njGlobal_) * niGlobal_);
  lonV.resize(static_cast<size_t>(njGlobal_) * niGlobal_);
  latV.resize(static_cast<size_t>(njGlobal_) * niGlobal_);

  for (int jG = 0; jG < njGlobal_; ++jG) {
    const int jSuper = 2 * jG + 1;
    for (int iG = 0; iG < niGlobal_; ++iG) {
      const int n = jG * niGlobal_ + iG;
      lonU[n] = xFull[ jSuper      * nxp + (2 * iG + 2)];
      latU[n] = yFull[ jSuper      * nxp + (2 * iG + 2)];
      lonV[n] = xFull[(jSuper + 1) * nxp + (2 * iG + 1)];
      latV[n] = yFull[(jSuper + 1) * nxp + (2 * iG + 1)];
    }
  }
}

// -------------------------------------------------------------------------------------------------
// Read depth(ny, nx) and wet(ny, nx) from ocean_topog.nc — full global arrays,
// size [njGlobal * niGlobal].  `wet` is MOM6's authoritative land/sea mask (1=ocean, 0=land).
// If `wet` is absent (e.g., synthetic test grids), falls back to depth > minimumDepth_.
void GeometryMOM6::readTopog(const std::string & path,
                              std::vector<double> & depth,
                              std::vector<double> & wet) const
{
  int ncid;
  if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
    throw eckit::CantOpenFile(path, Here());

  const size_t n = static_cast<size_t>(njGlobal_) * niGlobal_;
  depth.resize(n);
  wet.resize(n);
  int vid;
  nc_inq_varid(ncid, "depth", &vid); nc_get_var_double(ncid, vid, depth.data());
  // wet is optional; fall back to depth > minimumDepth_ when not present
  if (nc_inq_varid(ncid, "wet", &vid) == NC_NOERR) {
    nc_get_var_double(ncid, vid, wet.data());
  } else {
    for (size_t k = 0; k < n; ++k)
      wet[k] = (depth[k] > minimumDepth_) ? 1.0 : 0.0;
  }
  nc_close(ncid);
}

// -------------------------------------------------------------------------------------------------
// Build Atlas NodeColumns for the local MOM6 compute domain → mom6FunctionSpace_.
// Global arrays are passed in; local slice is extracted here.
void GeometryMOM6::buildMom6FunctionSpace(const eckit::mpi::Comm & comm,
                                           const std::vector<double> & lonGlobal,
                                           const std::vector<double> & latGlobal)
{
  using atlas::gidx_t;
  using atlas::idx_t;

  const int rank   = static_cast<int>(comm.rank());
  const int nOwned = iCount_ * jCount_;

  std::vector<double>  lons(nOwned), lats(nOwned);
  std::vector<int>     ghosts(nOwned, 0), partitions(nOwned, rank);
  std::vector<gidx_t>  globalIdx(nOwned);
  std::vector<idx_t>   remoteIdx(nOwned);

  for (int jL = 0; jL < jCount_; ++jL) {
    for (int iL = 0; iL < iCount_; ++iL) {
      const int n  = jL * iCount_ + iL;
      const int iG = iStart_ - 1 + iL;
      const int jG = jStart_ - 1 + jL;
      lons[n]      = lonGlobal[jG * niGlobal_ + iG];
      lats[n]      = latGlobal[jG * niGlobal_ + iG];
      globalIdx[n] = static_cast<gidx_t>(jG * niGlobal_ + iG + 1);  // 1-based
      remoteIdx[n] = static_cast<idx_t>(n + 1);                      // 1-based
    }
  }

  eckit::LocalConfiguration meshConf;
  meshConf.set("mpi_comm", comm.name());

  const atlas::mesh::MeshBuilder meshBuilder;
  atlas::Mesh mesh = meshBuilder(
      lons, lats, ghosts,
      globalIdx, remoteIdx, /*remote_index_base=*/1, partitions,
      /*tri_boundary_nodes=*/{}, /*tri_global_indices=*/{},
      /*quad_boundary_nodes=*/{}, /*quad_global_indices=*/{},
      meshConf);

  mom6FunctionSpace_ = atlas::functionspace::NodeColumns(
      mesh, atlas::util::Config("mpi_comm", comm.name()));
}

// -------------------------------------------------------------------------------------------------
// Build Atlas NodeColumns for the JEDI unstructured partition → functionSpace_ (base class).
//
// Algorithm:
//   1. All ranks build the same global active-point list (ocean + land fringe) — no MPI needed,
//      deterministic from the global mask.
//   2. Atlas "equal_regions" partitioner assigns each active point to a rank based on its
//      lon/lat, giving geographically compact (latitude-band) ownership.
//   3. Ghost nodes: active-point neighbours of owned nodes that belong to other JEDI ranks.
//   4. jediPoints_ records (iG, jG) for each node (owned first, then ghost).
//   5. Atlas MeshBuilder → NodeColumns.
//
// Periodicity: i wraps (east-west); j is clamped (no tripolar fold for now).
void GeometryMOM6::buildJediFunctionSpace(const eckit::mpi::Comm & comm,
                                           const std::vector<double> & wetGlobal,
                                           const std::vector<double> & lonGlobal,
                                           const std::vector<double> & latGlobal)
{
  using atlas::gidx_t;
  using atlas::idx_t;

  // Periodic ocean test using MOM6's authoritative wet mask (i wraps, j clamped)
  auto isOcean = [&](int i, int j) -> bool {
    i = (i % niGlobal_ + niGlobal_) % niGlobal_;
    if (j < 0 || j >= njGlobal_) return false;
    return wetGlobal[j * niGlobal_ + i] > 0.5;
  };

  // --- Step 1: global active-point list (row-major scan order) ---
  struct APoint { int iG, jG; };
  std::vector<APoint> active;
  active.reserve(niGlobal_ * njGlobal_);

  for (int jG = 0; jG < njGlobal_; ++jG) {
    for (int iG = 0; iG < niGlobal_; ++iG) {
      const bool ocean  = isOcean(iG, jG);
      const bool fringe = !ocean && (isOcean(iG - 1, jG) || isOcean(iG + 1, jG) ||
                                     isOcean(iG, jG - 1) || isOcean(iG, jG + 1));
      if (ocean || fringe)
        active.push_back({iG, jG});
    }
  }
  nActiveGlobal_ = static_cast<int>(active.size());

  // --- Step 2: Atlas geographic partition ---
  // Build an UnstructuredGrid from the active-point coordinates and let the
  // "equal_regions" partitioner assign geographically compact latitude-band ownership.
  std::vector<atlas::PointXY> pts(nActiveGlobal_);
  for (int k = 0; k < nActiveGlobal_; ++k) {
    const int gIdx = active[k].jG * niGlobal_ + active[k].iG;
    pts[k] = atlas::PointXY(lonGlobal[gIdx], latGlobal[gIdx]);
  }
  const atlas::UnstructuredGrid ugrid(pts);

  const int npes = static_cast<int>(comm.size());
  const int rank = static_cast<int>(comm.rank());

  const atlas::grid::Distribution dist =
      atlas::grid::Partitioner("equal_regions", npes).partition(ugrid);

  // partOf[k] = JEDI rank owning active point k
  std::vector<int> partOf(nActiveGlobal_);
  for (int k = 0; k < nActiveGlobal_; ++k)
    partOf[k] = dist.partition(k);

  ownedCount_ = static_cast<int>(dist.nb_pts()[rank]);
  ownedMin_   = static_cast<int>(dist.min_pts());
  ownedMax_   = static_cast<int>(dist.max_pts());

  // Fast lookup: flat global index (jG*NI+iG) → index in active[]
  std::unordered_map<int, int> activeMap;
  activeMap.reserve(nActiveGlobal_);
  for (int k = 0; k < nActiveGlobal_; ++k)
    activeMap[active[k].jG * niGlobal_ + active[k].iG] = k;

  // --- Step 3: identify ghost nodes ---
  // First collect owned points, then scan their neighbours for cross-rank active points.
  std::vector<int> ghostIdxVec;
  std::unordered_map<int, int> ghostSeen;   // active-list index → position in ghostIdxVec

  for (int k = 0; k < nActiveGlobal_; ++k) {
    if (partOf[k] != rank) continue;
    const int iG = active[k].iG;
    const int jG = active[k].jG;
    const int nbI[4] = { (iG - 1 + niGlobal_) % niGlobal_,
                         (iG + 1) % niGlobal_, iG, iG };
    const int nbJ[4] = { jG, jG, jG - 1, jG + 1 };
    for (int d = 0; d < 4; ++d) {
      if (nbJ[d] < 0 || nbJ[d] >= njGlobal_) continue;
      auto it = activeMap.find(nbJ[d] * niGlobal_ + nbI[d]);
      if (it == activeMap.end()) continue;
      const int neighborK = it->second;
      if (partOf[neighborK] == rank) continue;
      if (ghostSeen.emplace(neighborK, static_cast<int>(ghostIdxVec.size())).second)
        ghostIdxVec.push_back(neighborK);
    }
  }

  // --- Step 4: populate jediPoints_ (owned first, then ghost) ---
  const int nGhost = static_cast<int>(ghostIdxVec.size());
  const int nTotal = ownedCount_ + nGhost;
  jediPoints_.clear();
  jediPoints_.reserve(nTotal);
  for (int k = 0; k < nActiveGlobal_; ++k)
    if (partOf[k] == rank)
      jediPoints_.push_back({active[k].iG, active[k].jG});
  for (int gk : ghostIdxVec)
    jediPoints_.push_back({active[gk].iG, active[gk].jG});

  // --- Step 5: build Atlas mesh ---
  std::vector<double>  lons(nTotal), lats(nTotal);
  std::vector<int>     ghosts(nTotal, 0), partitions(nTotal, rank);
  std::vector<gidx_t>  globalIdx(nTotal);
  std::vector<idx_t>   remoteIdx(nTotal);

  for (int n = 0; n < nTotal; ++n) {
    const int iG = jediPoints_[n].first;
    const int jG = jediPoints_[n].second;
    lons[n]      = lonGlobal[jG * niGlobal_ + iG];
    lats[n]      = latGlobal[jG * niGlobal_ + iG];
    globalIdx[n] = static_cast<gidx_t>(jG * niGlobal_ + iG + 1);   // 1-based
    remoteIdx[n] = static_cast<idx_t>(n + 1);                       // 1-based
  }
  for (int g = 0; g < nGhost; ++g) {
    const int n = ownedCount_ + g;
    ghosts[n]     = 1;
    partitions[n] = partOf[ghostIdxVec[g]];
  }

  eckit::LocalConfiguration meshConf;
  meshConf.set("mpi_comm", comm.name());

  const atlas::mesh::MeshBuilder meshBuilder;
  atlas::Mesh mesh = meshBuilder(
      lons, lats, ghosts,
      globalIdx, remoteIdx, /*remote_index_base=*/1, partitions,
      /*tri_boundary_nodes=*/{}, /*tri_global_indices=*/{},
      /*quad_boundary_nodes=*/{}, /*quad_global_indices=*/{},
      meshConf);

  functionSpace_ = atlas::functionspace::NodeColumns(
      mesh, atlas::util::Config("mpi_comm", comm.name()));

  oops::Log::info() << "GeometryMOM6 JEDI space: rank " << rank
                    << " nActiveGlobal=" << nActiveGlobal_
                    << " ownedCount="    << ownedCount_
                    << " ghostCount="    << nGhost << std::endl;
}

// -------------------------------------------------------------------------------------------------
// Populate mom6Fields_ on mom6FunctionSpace_ with geometry data for the local MOM6 compute domain.
// All input arrays are global (size njGlobal*niGlobal, row-major); local slice is extracted here.
void GeometryMOM6::buildMom6Fields(const std::vector<double> & lonGlobal,
                                    const std::vector<double> & latGlobal,
                                    const std::vector<double> & depthGlobal,
                                    const std::vector<double> & wetGlobal,
                                    const std::vector<double> & dxTGlobal,
                                    const std::vector<double> & dyTGlobal,
                                    const std::vector<double> & areaTGlobal,
                                    const std::vector<double> & lonUGlobal,
                                    const std::vector<double> & latUGlobal,
                                    const std::vector<double> & lonVGlobal,
                                    const std::vector<double> & latVGlobal)
{
  const int npts = iCount_ * jCount_;
  mom6Fields_ = atlas::FieldSet();

  auto addField = [&](const std::string & name) -> atlas::Field {
    atlas::Field f = mom6FunctionSpace_.createField<double>(
        atlas::option::name(name) | atlas::option::levels(1));
    mom6Fields_.add(f);
    return f;
  };

  atlas::Field fOwned = mom6FunctionSpace_.createField<int>(
      atlas::option::name("owned") | atlas::option::levels(1));
  auto vOwned = atlas::array::make_view<int, 2>(fOwned);
  for (int n = 0; n < npts; ++n) vOwned(n, 0) = 1;
  mom6Fields_.add(fOwned);

  atlas::Field fLon   = addField("lon");
  atlas::Field fLat   = addField("lat");
  atlas::Field fMask  = addField("mask2d");
  atlas::Field fDepth = addField("depth");
  atlas::Field fDxT   = addField("dxT");
  atlas::Field fDyT   = addField("dyT");
  atlas::Field fAreaT = addField("areaT");
  atlas::Field fLonU  = addField("lonu");
  atlas::Field fLatU  = addField("latu");
  atlas::Field fLonV  = addField("lonv");
  atlas::Field fLatV  = addField("latv");

  auto vLon   = atlas::array::make_view<double, 2>(fLon);
  auto vLat   = atlas::array::make_view<double, 2>(fLat);
  auto vMask  = atlas::array::make_view<double, 2>(fMask);
  auto vDepth = atlas::array::make_view<double, 2>(fDepth);
  auto vDxT   = atlas::array::make_view<double, 2>(fDxT);
  auto vDyT   = atlas::array::make_view<double, 2>(fDyT);
  auto vAreaT = atlas::array::make_view<double, 2>(fAreaT);
  auto vLonU  = atlas::array::make_view<double, 2>(fLonU);
  auto vLatU  = atlas::array::make_view<double, 2>(fLatU);
  auto vLonV  = atlas::array::make_view<double, 2>(fLonV);
  auto vLatV  = atlas::array::make_view<double, 2>(fLatV);

  for (int jL = 0; jL < jCount_; ++jL) {
    for (int iL = 0; iL < iCount_; ++iL) {
      const int n    = jL * iCount_ + iL;
      const int gIdx = (jStart_ - 1 + jL) * niGlobal_ + (iStart_ - 1 + iL);
      vLon(n, 0)   = lonGlobal[gIdx];
      vLat(n, 0)   = latGlobal[gIdx];
      vDepth(n, 0) = depthGlobal[gIdx];
      vMask(n, 0)  = wetGlobal[gIdx];
      vDxT(n, 0)   = dxTGlobal[gIdx];
      vDyT(n, 0)   = dyTGlobal[gIdx];
      vAreaT(n, 0) = areaTGlobal[gIdx];
      vLonU(n, 0)  = lonUGlobal[gIdx];
      vLatU(n, 0)  = latUGlobal[gIdx];
      vLonV(n, 0)  = lonVGlobal[gIdx];
      vLatV(n, 0)  = latVGlobal[gIdx];
    }
  }
}

// -------------------------------------------------------------------------------------------------
// Populate fields_ on the JEDI unstructured functionSpace_ from global arrays.
// Node ordering follows jediPoints_ (owned nodes first, then ghost nodes).
void GeometryMOM6::buildFields(const std::vector<double> & lonGlobal,
                                const std::vector<double> & latGlobal,
                                const std::vector<double> & depthGlobal,
                                const std::vector<double> & wetGlobal,
                                const std::vector<double> & dxTGlobal,
                                const std::vector<double> & dyTGlobal,
                                const std::vector<double> & areaTGlobal,
                                const std::vector<double> & lonUGlobal,
                                const std::vector<double> & latUGlobal,
                                const std::vector<double> & lonVGlobal,
                                const std::vector<double> & latVGlobal)
{
  const int npts = static_cast<int>(jediPoints_.size());  // owned + ghost
  fields_ = atlas::FieldSet();

  auto addField = [&](const std::string & name) -> atlas::Field {
    atlas::Field f = functionSpace_.createField<double>(
        atlas::option::name(name) | atlas::option::levels(1));
    fields_.add(f);
    return f;
  };

  // "owned": 1 for owned nodes, 0 for ghost nodes
  atlas::Field fOwned = functionSpace_.createField<int>(
      atlas::option::name("owned") | atlas::option::levels(1));
  auto vOwned = atlas::array::make_view<int, 2>(fOwned);
  for (int n = 0; n < npts; ++n) vOwned(n, 0) = (n < ownedCount_) ? 1 : 0;
  fields_.add(fOwned);

  atlas::Field fLon   = addField("lon");
  atlas::Field fLat   = addField("lat");
  atlas::Field fMask  = addField("mask2d");
  atlas::Field fDepth = addField("depth");
  atlas::Field fDxT   = addField("dxT");
  atlas::Field fDyT   = addField("dyT");
  atlas::Field fAreaT = addField("areaT");
  atlas::Field fLonU  = addField("lonu");
  atlas::Field fLatU  = addField("latu");
  atlas::Field fLonV  = addField("lonv");
  atlas::Field fLatV  = addField("latv");

  auto vLon   = atlas::array::make_view<double, 2>(fLon);
  auto vLat   = atlas::array::make_view<double, 2>(fLat);
  auto vMask  = atlas::array::make_view<double, 2>(fMask);
  auto vDepth = atlas::array::make_view<double, 2>(fDepth);
  auto vDxT   = atlas::array::make_view<double, 2>(fDxT);
  auto vDyT   = atlas::array::make_view<double, 2>(fDyT);
  auto vAreaT = atlas::array::make_view<double, 2>(fAreaT);
  auto vLonU  = atlas::array::make_view<double, 2>(fLonU);
  auto vLatU  = atlas::array::make_view<double, 2>(fLatU);
  auto vLonV  = atlas::array::make_view<double, 2>(fLonV);
  auto vLatV  = atlas::array::make_view<double, 2>(fLatV);

  for (int n = 0; n < npts; ++n) {
    const int iG   = jediPoints_[n].first;
    const int jG   = jediPoints_[n].second;
    const int gIdx = jG * niGlobal_ + iG;
    vLon(n, 0)   = lonGlobal[gIdx];
    vLat(n, 0)   = latGlobal[gIdx];
    vDepth(n, 0) = depthGlobal[gIdx];
    vMask(n, 0)  = wetGlobal[gIdx];
    vDxT(n, 0)   = dxTGlobal[gIdx];
    vDyT(n, 0)   = dyTGlobal[gIdx];
    vAreaT(n, 0) = areaTGlobal[gIdx];
    vLonU(n, 0)  = lonUGlobal[gIdx];
    vLatU(n, 0)  = latUGlobal[gIdx];
    vLonV(n, 0)  = lonVGlobal[gIdx];
    vLatV(n, 0)  = latVGlobal[gIdx];
  }
}

// -------------------------------------------------------------------------------------------------
// Build scatter/gather map: for each JEDI-owned point, record the owning MOM6 rank and
// the row-major local index on that rank.  Deterministic — no MPI needed.
void GeometryMOM6::buildScatterMap()
{
  scatterMap_.mom6Rank.resize(ownedCount_);
  scatterMap_.mom6LocalIdx.resize(ownedCount_);

  for (int n = 0; n < ownedCount_; ++n) {
    const int iG = jediPoints_[n].first;
    const int jG = jediPoints_[n].second;

    // Find MOM6 piX (i-rank) and piY (j-rank) owning (iG, jG)
    const int baseX  = niGlobal_ / layoutX_;
    const int extraX = niGlobal_ % layoutX_;
    const int piX    = (iG < extraX * (baseX + 1))
                       ? iG / (baseX + 1)
                       : extraX + (iG - extraX * (baseX + 1)) / baseX;

    const int baseY  = njGlobal_ / layoutY_;
    const int extraY = njGlobal_ % layoutY_;
    const int piY    = (jG < extraY * (baseY + 1))
                       ? jG / (baseY + 1)
                       : extraY + (jG - extraY * (baseY + 1)) / baseY;

    scatterMap_.mom6Rank[n] = piX * layoutY_ + piY;

    // 0-based start of that tile
    int iStart0 = 0;
    for (int k = 0; k < piX; ++k) iStart0 += computeExtent(niGlobal_, layoutX_, k);
    int jStart0 = 0;
    for (int k = 0; k < piY; ++k) jStart0 += computeExtent(njGlobal_, layoutY_, k);

    scatterMap_.mom6LocalIdx[n] = (jG - jStart0) * computeExtent(niGlobal_, layoutX_, piX)
                                + (iG - iStart0);
  }
}

// -------------------------------------------------------------------------------------------------
GeometryMOM6::GeometryMOM6(const eckit::Configuration & conf,
                           const eckit::mpi::Comm & comm)
{
  oops::Log::trace() << "GeometryMOM6 constructor starting" << std::endl;

  const std::string inputDir = conf.getString("input_dir", ".");

  // 1. Read grid parameters from the MOM_input sub-configuration
  const eckit::LocalConfiguration momConf(conf, "MOM_input");
  niGlobal_     = momConf.getInt("NIGLOBAL");
  njGlobal_     = momConf.getInt("NJGLOBAL");
  numLevels_    = momConf.getInt("NZ");
  minimumDepth_ = momConf.getDouble("MINIMUM_DEPTH");
  const std::vector<int> layout = momConf.getIntVector("LAYOUT");
  layoutX_ = layout[0];
  layoutY_ = layout[1];

  oops::Log::debug() << "GeometryMOM6: NI=" << niGlobal_ << " NJ=" << njGlobal_
                     << " NZ=" << numLevels_
                     << " layout=(" << layoutX_ << "," << layoutY_ << ")" << std::endl;

  // 2. Compute local MOM6 domain extent for this rank
  //    MOM6 convention: rank = piX * layoutY_ + piY
  const int rank = static_cast<int>(comm.rank());
  const int piX  = rank / layoutY_;
  const int piY  = rank % layoutY_;

  iCount_ = computeExtent(niGlobal_, layoutX_, piX);
  jCount_ = computeExtent(njGlobal_, layoutY_, piY);

  iStart_ = 1;
  for (int k = 0; k < piX; ++k) iStart_ += computeExtent(niGlobal_, layoutX_, k);
  jStart_ = 1;
  for (int k = 0; k < piY; ++k) jStart_ += computeExtent(njGlobal_, layoutY_, k);

  oops::Log::debug() << "GeometryMOM6 rank " << rank
                     << ": iStart=" << iStart_ << " iCount=" << iCount_
                     << " jStart=" << jStart_ << " jCount=" << jCount_ << std::endl;

  // 3. Read global grid arrays (all ranks read identically; each is njGlobal*niGlobal)
  const std::string hgridPath = inputDir + "/INPUT/ocean_hgrid.nc";

  std::vector<double> lon, lat;
  readHgridLonLat(hgridPath, lon, lat);

  std::vector<double> dxT, dyT, areaT;
  readHgridMetrics(hgridPath, dxT, dyT, areaT);

  std::vector<double> lonU, latU, lonV, latV;
  readHgridUVLonLat(hgridPath, lonU, latU, lonV, latV);

  std::vector<double> depth, wet;
  readTopog(inputDir + "/INPUT/ocean_topog.nc", depth, wet);

  // 4. MOM6 structured function space and fields (local compute domain)
  buildMom6FunctionSpace(comm, lon, lat);
  buildMom6Fields(lon, lat, depth, wet, dxT, dyT, areaT, lonU, latU, lonV, latV);

  // 5. JEDI unstructured function space and fields (ocean + fringe, load-balanced)
  buildJediFunctionSpace(comm, wet, lon, lat);
  buildFields(lon, lat, depth, wet, dxT, dyT, areaT, lonU, latU, lonV, latV);

  // 6. Scatter/gather map: JEDI ↔ MOM6
  buildScatterMap();

  // 7. Optionally save grids to NetCDF / debug files
  if (conf.has("save grid to"))
    saveGrid(conf.getString("save grid to"), comm);
  if (conf.has("save structured grid to"))
    saveStructuredGrid(conf.getString("save structured grid to"), comm);
  if (conf.has("save debug mesh to"))
    saveDebugMesh(conf.getString("save debug mesh to"), comm);

  oops::Log::trace() << "GeometryMOM6 constructor done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
// Gather all MOM6-structured geometry fields from all ranks and write a (nj, ni) NetCDF
// file on rank 0.  Uses mom6Fields_ which holds the local MOM6 compute-domain data.
void GeometryMOM6::saveStructuredGrid(const std::string & filename,
                                       const eckit::mpi::Comm & comm) const
{
  const int npes = static_cast<int>(comm.size());
  const int root = 0;

  const std::vector<std::string> fieldNames =
      {"lon", "lat", "dxT", "dyT", "areaT", "depth", "mask2d",
       "lonu", "latu", "lonv", "latv"};

  // Step 1: gather per-rank local sizes on root
  const int localSize = jCount_ * iCount_;
  std::vector<int> recvcounts(npes, 0);
  comm.gather(localSize, recvcounts, root);

  std::vector<int> displs;
  if (comm.rank() == root) {
    displs.resize(npes);
    displs[0] = 0;
    for (int p = 1; p < npes; ++p)
      displs[p] = displs[p - 1] + recvcounts[p - 1];
  }

  // Step 2: on root, precompute (iStart, jStart, iCount, jCount) for every rank
  struct DomInfo { int iStart, jStart, iCount, jCount; };
  std::vector<DomInfo> allDoms;
  if (comm.rank() == root) {
    allDoms.resize(npes);
    for (int p = 0; p < npes; ++p) {
      const int px = p / layoutY_;
      const int py = p % layoutY_;
      allDoms[p].iCount = computeExtent(niGlobal_, layoutX_, px);
      allDoms[p].jCount = computeExtent(njGlobal_, layoutY_, py);
      allDoms[p].iStart = 1;
      for (int k = 0; k < px; ++k)
        allDoms[p].iStart += computeExtent(niGlobal_, layoutX_, k);
      allDoms[p].jStart = 1;
      for (int k = 0; k < py; ++k)
        allDoms[p].jStart += computeExtent(njGlobal_, layoutY_, k);
    }
  }

  // Step 3: create NetCDF file on root
  int ncid = -1;
  std::vector<int> varids(fieldNames.size(), -1);
  if (comm.rank() == root) {
    if (nc_create(filename.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid) != NC_NOERR)
      throw eckit::CantOpenFile(filename, Here());
    int nj_dim, ni_dim;
    nc_def_dim(ncid, "nj", static_cast<size_t>(njGlobal_), &nj_dim);
    nc_def_dim(ncid, "ni", static_cast<size_t>(niGlobal_), &ni_dim);
    const int dims[2] = {nj_dim, ni_dim};
    for (size_t f = 0; f < fieldNames.size(); ++f)
      nc_def_var(ncid, fieldNames[f].c_str(), NC_DOUBLE, 2, dims, &varids[f]);
    nc_enddef(ncid);
  }

  // Step 4: for each field, gather and write
  std::vector<double> recvBuf;
  if (comm.rank() == root) {
    const int total = displs.back() + recvcounts.back();
    recvBuf.resize(total);
  }

  for (size_t f = 0; f < fieldNames.size(); ++f) {
    auto view = atlas::array::make_view<double, 2>(mom6Fields_.field(fieldNames[f]));
    std::vector<double> localData(localSize);
    for (int n = 0; n < localSize; ++n) localData[n] = view(n, 0);

    comm.gatherv(localData, recvBuf, recvcounts, displs, root);

    if (comm.rank() == root) {
      std::vector<double> global(static_cast<size_t>(njGlobal_) * niGlobal_, 0.0);
      for (int p = 0; p < npes; ++p) {
        const DomInfo & d = allDoms[p];
        const int off = displs[p];
        for (int jL = 0; jL < d.jCount; ++jL)
          for (int iL = 0; iL < d.iCount; ++iL) {
            const int jG = d.jStart - 1 + jL;
            const int iG = d.iStart - 1 + iL;
            global[jG * niGlobal_ + iG] = recvBuf[off + jL * d.iCount + iL];
          }
      }
      nc_put_var_double(ncid, varids[f], global.data());
    }
  }

  if (comm.rank() == root) nc_close(ncid);
}

// -------------------------------------------------------------------------------------------------
// Save the JEDI geometry fields to NetCDF using oops Atlas IO helpers.
void GeometryMOM6::saveGrid(const std::string & filename,
                             const eckit::mpi::Comm & comm) const
{
  atlas::FieldSet toWrite;
  for (const char * name : {"depth", "mask2d", "dxT", "dyT", "areaT",
                            "lonu", "latu", "lonv", "latv"})
    toWrite.add(fields_.field(name));

  std::string filepath = filename;
  const std::string ext = ".nc";
  if (filepath.size() > ext.size() &&
      filepath.compare(filepath.size() - ext.size(), ext.size(), ext) == 0)
    filepath.erase(filepath.size() - ext.size());

  eckit::LocalConfiguration conf;
  conf.set("filepath", filepath);
  util::writeFieldSet(comm, conf, toWrite);
}

// -------------------------------------------------------------------------------------------------
// Write debug output for inspecting Atlas domain decomposition and halo structure.
void GeometryMOM6::saveDebugMesh(const std::string & prefix,
                                  const eckit::mpi::Comm & comm) const
{
  const int rank = static_cast<int>(comm.rank());
  const atlas::functionspace::NodeColumns fspace(functionSpace_);
  const atlas::Mesh & mesh = fspace.mesh();

  {
    atlas::output::Gmsh gmsh(prefix + ".msh",
        atlas::util::Config("ghost",       true)
      | atlas::util::Config("info",        true)
      | atlas::util::Config("coordinates", "lonlat"));
    gmsh.write(mesh);
    gmsh.write(fields_.field("mask2d"), functionSpace_);
  }

  {
    const std::string ncFile =
        prefix + "_rank" + std::to_string(rank) + ".nc";

    const atlas::mesh::Nodes & nodes = mesh.nodes();
    const int nNodes = nodes.size();

    auto lonlatView = atlas::array::make_view<double,        2>(nodes.lonlat());
    auto ghostView  = atlas::array::make_view<int,           1>(nodes.ghost());
    auto partView   = atlas::array::make_view<int,           1>(nodes.partition());
    auto gidxView   = atlas::array::make_view<atlas::gidx_t, 1>(nodes.global_index());
    auto maskView   = atlas::array::make_view<double,        2>(fields_.field("mask2d"));

    std::vector<double>    lons(nNodes), lats(nNodes), mask(nNodes);
    std::vector<int>       ghost(nNodes), part(nNodes);
    std::vector<long long> gidx(nNodes);
    for (int n = 0; n < nNodes; ++n) {
      lons[n]  = lonlatView(n, 0);
      lats[n]  = lonlatView(n, 1);
      ghost[n] = ghostView(n);
      part[n]  = partView(n);
      gidx[n]  = static_cast<long long>(gidxView(n));
      mask[n]  = maskView(n, 0);
    }

    int ncid;
    if (nc_create(ncFile.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid) != NC_NOERR)
      throw eckit::CantOpenFile(ncFile, Here());

    int node_dim;
    nc_def_dim(ncid, "node", static_cast<size_t>(nNodes), &node_dim);

    int vid_lon, vid_lat, vid_ghost, vid_part, vid_gidx, vid_rank, vid_mask;
    nc_def_var(ncid, "lon",          NC_DOUBLE, 1, &node_dim, &vid_lon);
    nc_def_var(ncid, "lat",          NC_DOUBLE, 1, &node_dim, &vid_lat);
    nc_def_var(ncid, "ghost",        NC_INT,    1, &node_dim, &vid_ghost);
    nc_def_var(ncid, "partition",    NC_INT,    1, &node_dim, &vid_part);
    nc_def_var(ncid, "global_index", NC_INT64,  1, &node_dim, &vid_gidx);
    nc_def_var(ncid, "mask2d",       NC_DOUBLE, 1, &node_dim, &vid_mask);
    nc_def_var(ncid, "rank",         NC_INT,    0, nullptr,   &vid_rank);
    nc_enddef(ncid);

    nc_put_var_double  (ncid, vid_lon,   lons.data());
    nc_put_var_double  (ncid, vid_lat,   lats.data());
    nc_put_var_int     (ncid, vid_ghost, ghost.data());
    nc_put_var_int     (ncid, vid_part,  part.data());
    nc_put_var_longlong(ncid, vid_gidx,  gidx.data());
    nc_put_var_double  (ncid, vid_mask,  mask.data());
    nc_put_var_int     (ncid, vid_rank,  &rank);
    nc_close(ncid);

    oops::Log::info() << "GeometryMOM6::saveDebugMesh: rank " << rank
                      << " -> " << ncFile << std::endl;
  }
}

// -------------------------------------------------------------------------------------------------
void GeometryMOM6::print(std::ostream & os) const
{
  const int    npes     = layoutX_ * layoutY_;
  const double avg      = static_cast<double>(nActiveGlobal_) / npes;
  const double imbalPct = (ownedMax_ > 0)
                          ? 100.0 * (ownedMax_ - ownedMin_) / avg
                          : 0.0;
  os << "\n"
     << "  +----- MOM6 Geometry ------------------------------------------+\n"
     << "  |  Grid      : " << niGlobal_ << " x " << njGlobal_ << " x " << numLevels_
                            << "  (NI x NJ x NZ)\n"
     << "  |  MOM6 layout : " << layoutX_ << " x " << layoutY_
                              << "  (" << npes << " MPI tasks)\n"
     << "  |  Min depth   : " << minimumDepth_ << " m\n"
     << "  |  ---- JEDI unstructured partition (equal_regions) ----------\n"
     << "  |  Active pts  : " << nActiveGlobal_
                              << "  (ocean + land fringe)\n"
     << "  |  Per rank    : avg " << static_cast<int>(avg)
                                  << "  min " << ownedMin_
                                  << "  max " << ownedMax_
                                  << "  imbalance "
                                  << std::fixed << std::setprecision(1)
                                  << imbalPct << "%\n"
     << "  |  This rank   : " << ownedCount_ << " owned pts\n"
     << "  +--------------------------------------------------------------+\n";
}

// -------------------------------------------------------------------------------------------------
eckit::LocalConfiguration GeometryMOM6::gridSpecific() const
{
  eckit::LocalConfiguration conf;
  conf.set("grid_type",      "mom6");
  conf.set("ni_global",      niGlobal_);
  conf.set("nj_global",      njGlobal_);
  conf.set("nz",             numLevels_);
  conf.set("layout_x",       layoutX_);
  conf.set("layout_y",       layoutY_);
  conf.set("minimum_depth",  minimumDepth_);
  conf.set("n_active_global", nActiveGlobal_);
  return conf;
}

// -------------------------------------------------------------------------------------------------
} // namespace ijedi

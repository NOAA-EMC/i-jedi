/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

 // A geometry class that only needs lat/lon/elevation information--for land (snow/soil) DA over independent grids.
 // Written with the intent of being a minimal geometry class that can be used for land DA over independent grids, without requiring a full model geometry.
 // Written with the help of Gemeni, fully reviewed and tested by land DA team.
  
#include <netcdf.h>

#include <string>
#include <vector>
#include <iomanip>  // Required for std::setprecision and std::fixed
#include <iostream>
#include <utility>

#include "atlas/array.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/mesh/actions/BuildHalo.h"

#include "atlas/meshgenerator.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
//#include "atlas/mesh/actions/BuildDelaunay.h"
//#include "atlas/mesh/actions/Dispatch.h"

#include "atlas/output/Gmsh.h"
#include "atlas/util/Earth.h"
#include "atlas/util/KDTree.h"

#include "atlas/util/Config.h"
#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/util/Metadata.h"
#include "oops/util/FieldSetHelpers.h"
#include "atlas/util/Earth.h"
#include "atlas/util/KDTree.h"

#include "ijedi/Geometry/landVector/GeometryLandVector.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include "atlas/util/Point.h"

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <array>
#include "atlas/array.h"
#include "atlas/util/Constants.h"


#include <tuple>
#include "eckit/mpi/Comm.h"


static inline void nc_rc(const int return_code, const std::string & operation) {
  if (return_code != NC_NOERR) {
    std::string errMsg = "GeometryLandVector netCDF operation '" + operation +
                         "' failed with error: " + nc_strerror(return_code);
    oops::Log::error() << errMsg << std::endl;
    throw eckit::Exception(errMsg, Here());
  }
}

namespace ijedi {

// -----------------------------------------------------------------------------
// POD Struct for transferring halo node metadata between MPI ranks
// -----------------------------------------------------------------------------
struct GhostNodeData {
  atlas::gidx_t globalIdx;
  int tile;
  int i;
  int j;
  double lon;
  double lat;
  atlas::idx_t remoteIdx;
  int ownerRank;
};

// -----------------------------------------------------------------------------
// Helper: Pack (tile, i, j) coordinates into a 64-bit lookup hash key
// -----------------------------------------------------------------------------
inline uint64_t pack_tile_ij(int tile, int i, int j) {
  return (static_cast<uint64_t>(tile) << 32) |
         (static_cast<uint64_t>(i & 0xFFFF) << 16) |
         static_cast<uint64_t>(j & 0xFFFF);
}

struct Fv3Corner {
  int t;
  int i;
  int j;
};

struct Fv3StepState {
  int t;
  int i;
  int j;
  int di;
  int dj;
  bool valid;
};

/*TODO: move the quad nodes construction step here
// Construct quads across panels using FV3 seam neighbor resolution
void compute_quad_nodes(
    size_t ownedCount_,
    int N,
    const std::vector<double> lons,
    const std::vector<double> lats,
    const std::vector<int>& rank_tiles,
    const std::vector<int>& rank_i,
    const std::vector<int>& rank_j,
    const std::vector<atlas::gidx_t>& globalIdx,
    std::vector<std::array<atlas::gidx_t, 4>>& quadNodes,
    std::vector<atlas::gidx_t>& quadGidx,
    const eckit::mpi::Comm& comm)
{
  quadNodes.clear();
  quadGidx.clear();



}*/
 
//============================================================================================


// Initialize static member variables
std::map<std::string, std::shared_ptr<GeometryLandVector::SharedData>>
    GeometryLandVector::cache_;
std::mutex GeometryLandVector::mutex_;

//---------------cache clear usage e.g. in tests--------------
/*#include "gtest/gtest.h"
#include "ijedi/Geometry/landVector/GeometryLandVector.h"

class GeometryLandVectorTest : public ::testing::Test {
 protected:
  void TearDown() override {
    // Purge cached geometries between tests
    ijedi::GeometryLandVector::clearCache();
    EXPECT_EQ(ijedi::GeometryLandVector::cacheSize(), 0u);
  }
};*/

// -----------------------------------------------------------------------------
// Clear Cache Function: Call this in C++ Test Fixtures (e.g. TearDown)
// -----------------------------------------------------------------------------
void GeometryLandVector::clearCache() {
  std::lock_guard<std::mutex> lock(mutex_);
  // Clearing the map drops reference counts.
  // If no GeometryLandVector objects exist in the test harness,
  // underlying Atlas grid/PointCloud memory is immediately freed.
  cache_.clear();
}

// -----------------------------------------------------------------------------
// Cache Size Inspector (Useful for assertion checks in unit tests)
// -----------------------------------------------------------------------------
size_t GeometryLandVector::cacheSize() {
  std::lock_guard<std::mutex> lock(mutex_);
  return cache_.size();
}

// -----------------------------------------------------------------------------
// Cache Key Generator
// -----------------------------------------------------------------------------
std::string GeometryLandVector::makeCacheKey(const eckit::Configuration & config,
                                             const eckit::mpi::Comm & comm) {
  // Use unique key parameters driving the PointCloud (e.g. filename, resolution)
  std::string grid_resol = config.getString("grid_resol", "C384.mx025");
  return "landVectorGeom_" + grid_resol + "_" + comm.name() + "_" + std::to_string(comm.size());
}

// -----------------------------------------------------------------------------

// Primary Constructor
// -----------------------------------------------------------------------------
GeometryLandVector::GeometryLandVector(const eckit::Configuration &config,
                                       const eckit::mpi::Comm &comm,
                                       eckit::LocalConfiguration & geomVariables,
                                       atlas::FunctionSpace &functionSpace,
                                       atlas::FieldSet &geomFields,
                                       bool &levelsAreTopDown, int &numLevels)
{
  // Read toggle flag from YAML (defaults to true if key is omitted)
  const bool useCache = config.getBool("geometry_cache", true);

  if (useCache) {
    std::string key = makeCacheKey(config, comm);

    std::lock_guard<std::mutex> lock(mutex_);

    auto it = cache_.find(key);
    if (it != cache_.end()) {
      // 1. INSTANCE CREATED FROM SCRATCH MATCHES EXISTING CACHE:
      // Reuse existing shared data handle across distinct OOPS factory calls
      sharedData_ = it->second;
    } else {
      // 2. FIRST TIME THIS GEOMETRY SPECIFICATION IS CREATED:
      // Build expensive spatial structures once and store handle in cache map
      sharedData_ = std::make_shared<SharedData>(config, comm, geomVariables, functionSpace,
                                         geomFields, levelsAreTopDown, numLevels);
      cache_[key] = sharedData_;
    }
  } else {
    // 3. CACHE BYPASSED: Construct dedicated instance directly
    sharedData_ = std::make_shared<SharedData>(config, comm, geomVariables, functionSpace,
                                               geomFields, levelsAreTopDown, numLevels);
  }
  // 4: populate output parameters
  functionSpace = sharedData_->functionSpace_;
  geomFields = sharedData_->fields_;
  geomVariables = sharedData_->geomVariables_;
  levelsAreTopDown = sharedData_->levelsAreTopDown_;
  numLevels = sharedData_->numLevels_;
  //numberLevels_ = sharedData_->numLevels_;

}

// -----------------------------------------------------------------------------
GeometryLandVector::SharedData::SharedData(const eckit::Configuration &config,
                                       const eckit::mpi::Comm &comm,
                                       eckit::LocalConfiguration & /*geomVariables*/,
                                       atlas::FunctionSpace &functionSpace,
                                       atlas::FieldSet &geomFields,
                                       bool &levelsAreTopDown, int &numLevels)
    : comm_(comm) 
{
  oops::Log::trace() << "GeometryLandVector::GeometryLandVector starting..." << std::endl;
 
  numLevels_ = config.getInt("nlevels", 1);
  double R_halo = config.getDouble("halo_radius", 125.0); //km
  numLevels = numLevels_;
  levelsAreTopDown = config.getBool("levels_are_topdown", false);
  levelsAreTopDown_ = levelsAreTopDown;

  bool print_debug_info = config.getBool("print_debug_info", false);

  // 1. Read coordinates (lat, lon, elevation) from input  file
  std::string latlonFile = config.getString("latlon_file", "");
  std::string latName = config.getString("lat_name", "latitude");
  std::string lonName = config.getString("lon_name", "longitude");
  std::string elevName = config.getString("elev_name", "elevation");

  std::vector<double> lats;
  std::vector<double> lons; 
  std::vector<double> elevations;
  std::vector<int> global_tile;
  std::vector<int> global_i;
  std::vector<int> global_j;

  if(!latlonFile.empty() && config.has("lat_name") && config.has("lon_name") && config.has("elev_name")) {
    oops::Log::info() << "Reading lat/lon/elevation from file: " << latlonFile << std::endl;
    readLatLonElevFromFile(latlonFile, latName, lonName, elevName, lats, lons, elevations, global_tile, global_i, global_j);
  } else {
    throw eckit::BadParameter("GeometryLandVector: Must specify latlon_file in configuration.", Here());
  }

  const size_t my_rank = static_cast<size_t>(comm_.rank());
  const size_t nproc = static_cast<size_t>(comm_.size());
  oops::Log::info() << ":GeometryLandVector: rank  " << my_rank << " running with " << nproc << " procs" << std::endl;

  // 2. divide points among ranks evenly
  const size_t num_points = lats.size();
  std::vector<atlas::PointXY> global_pts(num_points);
  for (size_t i = 0; i < num_points; ++i) {
    global_pts[i] = atlas::PointXY(lons[i], lats[i]);
  }
  
  size_t points_per_rank = num_points / nproc;
  size_t remainder = num_points % nproc;
  
  size_t start_idx = my_rank * points_per_rank + std::min<size_t>(my_rank, remainder);
  size_t local_count = points_per_rank + (my_rank < remainder ? 1 : 0);

  std::vector<int> rank_owning(num_points, 0);
  // Assign balanced chunks across MPI ranks
  for (size_t i = 0; i < num_points; ++i) {
    rank_owning[i] = static_cast<int>((i * nproc) / num_points);
  }
  
  // 3. Assign rank specific tile/index info 
  std::vector<atlas::PointXY> rank_points;
  std::vector<int> rank_tiles;
  std::vector<int> rank_i;
  std::vector<int> rank_j;
  std::vector<double> rank_lons;
  std::vector<double> rank_lats;
  std::vector<int> ghosts;  //(ownedCount_,0); // 0 for owned, 1 for ghost
  std::vector<atlas::gidx_t> globalIdx; // Original 1-based global index
  std::vector<atlas::idx_t> remoteIdx;      // use 1-based
  std::vector<int> partitions;

  //atlas::idx_t rem_indx = 1;
  std::vector<int> ownerLocalIdx(num_points);  // each active point's local index on its owner rank.
  std::vector<int> ownerCount(nproc, 0);
  for (int i = 0; i < num_points; ++i){
    int owner_rank = rank_owning[i];  
    ownerLocalIdx[i] = ownerCount[owner_rank]++;
  }
  //atlas::idx_t rem_indx = 1;
  for (size_t g_idx = 0; g_idx < global_pts.size(); ++g_idx) {
    if (rank_owning[g_idx] == my_rank){  
        rank_points.push_back(global_pts[g_idx]);
        ghosts.push_back(0); // Owned
        globalIdx.push_back(g_idx + 1); 
	rank_lons.push_back(global_pts[g_idx].x());
        rank_lats.push_back(global_pts[g_idx].y());
	remoteIdx.push_back(static_cast<atlas::idx_t>(ownerLocalIdx[g_idx] + 1));  //rem_indx++); 

	rank_tiles.push_back(global_tile[g_idx]);
        rank_i.push_back(global_i[g_idx]);
        rank_j.push_back(global_j[g_idx]);
        partitions.push_back(my_rank);
    }
  }
    
  ownedCount_ = rank_points.size(); 
  if (print_debug_info) {
      std::vector<size_t> num_owned_vec(nproc);
      comm_.allGather(ownedCount_, num_owned_vec.begin(), num_owned_vec.end());
      auto min_no = *std::min_element(num_owned_vec.begin(), num_owned_vec.end());
      auto max_no = *std::max_element(num_owned_vec.begin(), num_owned_vec.end());
      oops::Log::debug() << "proc " << my_rank << " ownedCount_ " << ownedCount_ << " min/max ownedCount_ " 
      	  << min_no << " / " << max_no << std::endl;  
  }

  // 4. Build global spatial index to discover owners of missing boundary vertices
  std::unordered_map<uint64_t, size_t> globalCoordToIdx;
  globalCoordToIdx.reserve(num_points);
  for (size_t g = 0; g < num_points; ++g) {
    globalCoordToIdx[pack_tile_ij(global_tile[g], global_i[g], global_j[g])] = g;
  }

  // 4b. Precompute perimeter cells (i=1, i=N, j=1, j=N) and geometric seam lookup
  // (Used by both Step 5 for MPI ghost requests and Step 9 for quad generation)
  const int N = config.getInt("resolution", 384);
  const int idxMin = 1; //
  const int idxMax = (idxMin == 1) ? N : (N - 1);

  struct BoundaryCell {
    int t, i, j;
    double x, y, z;
  };
  std::vector<BoundaryCell> boundaryCells;
  for (size_t g = 0; g < num_points; ++g) {
    int bt = global_tile[g];
    int bi = global_i[g];
    int bj = global_j[g];
    if (bi == idxMin || bi == idxMax || bj == idxMin || bj == idxMax ) {
      double rlon = lons[g] * M_PI / 180.0;
      double rlat = lats[g] * M_PI / 180.0;
      boundaryCells.push_back({
          bt, bi, bj,
          std::cos(rlat) * std::cos(rlon),
          std::cos(rlat) * std::sin(rlon),
          std::sin(rlat)});
    }
  }

  // Great-circle distance check on the unit sphere (avoids polar convergence and 0/360 wrap issues)
  // use value that safely allows stretched cubed-sphere corners while rejecting trans-oceanic jumps.
  const double maxEdgeDeg = config.getDouble("maximum_edge_size_deg", 3.0);
  const double maxChordSq = std::pow(2.0 * std::sin((maxEdgeDeg * M_PI / 180.0) / 2.0), 2);

  auto find_nearest_foreign_cell = [&](int t, int i, int j) -> std::pair<Fv3Corner, bool> {
    auto it = globalCoordToIdx.find(pack_tile_ij(t, i, j));
    if (it == globalCoordToIdx.end()) return {{}, false};
    size_t g = it->second;
    double rlon = lons[g] * M_PI / 180.0;
    double rlat = lats[g] * M_PI / 180.0;
    double x0 = std::cos(rlat) * std::cos(rlon);
    double y0 = std::cos(rlat) * std::sin(rlon);
    double z0 = std::sin(rlat);

    double bestD2 = maxChordSq;
    Fv3Corner bestCorner{};
    bool found = false;

    for (const auto& bc : boundaryCells) {
      if (bc.t == t) continue;
      double dx = x0 - bc.x;
      double dy = y0 - bc.y;
      double dz = z0 - bc.z;
      double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 < bestD2) {
        bestD2 = d2;
        bestCorner = {bc.t, bc.i, bc.j};
        found = true;
      }
    }
    return {bestCorner, found};
  };

  auto get_fv3_quad_corners = [&](int t, int i, int j)
      -> std::tuple<Fv3Corner, Fv3Corner, Fv3Corner, bool> {
    // Interior cell (both i+1 and j+1 are <= idxMax)
    if (i < idxMax && j < idxMax) {
      return {Fv3Corner{t, i + 1, j},
              Fv3Corner{t, i + 1, j + 1},
              Fv3Corner{t, i, j + 1},
              true};
    }
    // East tile seam (i == idxMax, j < idxMax)
    if (i == idxMax && j < idxMax) {
      auto [c_SE, ok_SE] = find_nearest_foreign_cell(t, i, j);
      auto [c_NE, ok_NE] = find_nearest_foreign_cell(t, i, j + 1);
      if (!ok_SE || !ok_NE) return {{}, {}, {}, false};
      return {c_SE, c_NE, Fv3Corner{t, i, j + 1}, true};
    }
    // North tile seam (j == idxMax, i < idxMax)
    if (j == idxMax && i < idxMax) {
      auto [c_NW, ok_NW] = find_nearest_foreign_cell(t, i, j);
      auto [c_NE, ok_NE] = find_nearest_foreign_cell(t, i + 1, j);
      if (!ok_NW || !ok_NE) return {{}, {}, {}, false};
      return {Fv3Corner{t, i + 1, j}, c_NE, c_NW, true};
    }
    // Exact tile corner (i == idxMax && j == idxMax): 3-tile triangle singularity
    return {{}, {}, {}, false};
  };  

  // 5. Determine neighbor nodes required to close boundary quads across FV3 seams
  std::vector<std::vector<atlas::gidx_t>> requests(nproc);

  for (size_t n = 0; n < ownedCount_; ++n) {
    int t = rank_tiles[n];
    int i = rank_i[n];
    int j = rank_j[n];

    auto [c_SE, c_NE, c_NW, valid] = get_fv3_quad_corners(t, i, j);
    if (!valid) continue;

    const std::array<Fv3Corner, 3> corners = {c_SE, c_NE, c_NW};
    for (const auto& c : corners) {
      uint64_t key = pack_tile_ij(c.t, c.i, c.j);
      auto it = globalCoordToIdx.find(key);
      if (it != globalCoordToIdx.end()) {
        size_t neighbor_g_idx = it->second;
        int owner = rank_owning[neighbor_g_idx];

        if (owner != my_rank) {
          atlas::gidx_t req_gid = static_cast<atlas::gidx_t>(neighbor_g_idx + 1);
          requests[owner].push_back(req_gid);
        }
      }
    }
  }

  // Deduplicate request lists per destination rank
  for (int r = 0; r < nproc; ++r) {
    std::sort(requests[r].begin(), requests[r].end());
    requests[r].erase(std::unique(requests[r].begin(), requests[r].end()), requests[r].end());
  }

  // 6. Non-blocking MPI halo exchange (Phase 1: Exchange node ID requests)
  std::vector<int> send_req_counts(nproc, 0);
  for (int r = 0; r < nproc; ++r) {
    send_req_counts[r] = static_cast<int>(requests[r].size());
  }

  std::vector<int> recv_req_counts(nproc, 0);
  comm_.allToAll(send_req_counts, recv_req_counts);

  std::vector<std::vector<atlas::gidx_t>> incoming_requests(nproc);
  std::vector<eckit::mpi::Request> req_stage1;

  for (int r = 0; r < nproc; ++r) {
    if (r != my_rank && recv_req_counts[r] > 0) {
      incoming_requests[r].resize(recv_req_counts[r]);
      req_stage1.push_back(comm_.iReceive(incoming_requests[r].data(),
                                         recv_req_counts[r],
                                         r, /*tag=*/100));
    }
  }

  for (int r = 0; r < nproc; ++r) {
    if (r != my_rank && send_req_counts[r] > 0) {
      req_stage1.push_back(comm_.iSend(requests[r].data(),
                                      send_req_counts[r],
                                      r, /*tag=*/100));
    }
  }

  comm_.waitAll(req_stage1);

  // 7. Non-blocking MPI halo exchange (Phase 2: Pack and send ghost node payloads)
  std::vector<std::vector<GhostNodeData>> outgoing_ghosts(nproc);
  for (int r = 0; r < nproc; ++r) {
    if (r == my_rank) continue;

    outgoing_ghosts[r].reserve(recv_req_counts[r]);
    for (atlas::gidx_t gid : incoming_requests[r]) {
      size_t g_idx = static_cast<size_t>(gid - 1);
      GhostNodeData node{};  //{} zero-initialize
      node.globalIdx = gid;
      node.tile      = global_tile[g_idx];
      node.i         = global_i[g_idx];
      node.j         = global_j[g_idx];
      node.lon       = global_pts[g_idx].x();
      node.lat       = global_pts[g_idx].y();
      node.remoteIdx = static_cast<atlas::idx_t>(ownerLocalIdx[g_idx] + 1);
      node.ownerRank = my_rank;
      outgoing_ghosts[r].push_back(node);
    }
  }

  std::vector<std::vector<GhostNodeData>> incoming_ghosts(nproc);
  std::vector<eckit::mpi::Request> req_stage2;

  for (int r = 0; r < nproc; ++r) {
    if (r != my_rank && send_req_counts[r] > 0) {
      incoming_ghosts[r].resize(send_req_counts[r]);
      req_stage2.push_back(comm_.iReceive(
          reinterpret_cast<char*>(incoming_ghosts[r].data()),
          static_cast<int>(send_req_counts[r] * sizeof(GhostNodeData)),
          r, /*tag=*/200));
    }
  }

  for (int r = 0; r < nproc; ++r) {
    if (r != my_rank && recv_req_counts[r] > 0) {
      req_stage2.push_back(comm_.iSend(
          reinterpret_cast<const char*>(outgoing_ghosts[r].data()),
          static_cast<int>(outgoing_ghosts[r].size() * sizeof(GhostNodeData)),
          r, /*tag=*/200));
    }
  }

  comm_.waitAll(req_stage2);

  // 8. Append received ghost/halo nodes to local structures
  // Keep track of all global indices already registered locally (both owned and ghosts)
  std::unordered_set<atlas::gidx_t> registeredGidx(globalIdx.begin(), globalIdx.end());
  for (int r = 0; r < nproc; ++r) {
    for (const auto& g : incoming_ghosts[r]) {
      // If this ghost node was already added by another quad or incoming message, SKIP IT
      if (!registeredGidx.insert(g.globalIdx).second) {
        continue;
      }
      rank_points.push_back(atlas::PointXY(g.lon, g.lat));
      rank_lons.push_back(g.lon);
      rank_lats.push_back(g.lat);
      rank_tiles.push_back(g.tile);
      rank_i.push_back(g.i);
      rank_j.push_back(g.j);

      ghosts.push_back(1); // 1 = Ghost node
      globalIdx.push_back(g.globalIdx);
      remoteIdx.push_back(g.remoteIdx);
      partitions.push_back(g.ownerRank);
    }
  }
  
  if (print_debug_info) {
      size_t nghost = std::accumulate(ghosts.begin(), ghosts.end(), 0);
      std::vector<size_t> num_ghost_vec(nproc);
      comm_.allGather(nghost, num_ghost_vec.begin(), num_ghost_vec.end());
      auto min_ng = *std::min_element(num_ghost_vec.begin(), num_ghost_vec.end());
      auto max_ng = *std::max_element(num_ghost_vec.begin(), num_ghost_vec.end());
      oops::Log::debug() << "rank " << my_rank << "num ghost cells " << nghost << " min/max num ghosts "
              << min_ng << " / " << max_ng << std::endl;  
  }

  // 9. Build quads across local and halo nodes using rotated seam transformations
  //    and great-circle edge length filtering
  std::vector<std::array<atlas::gidx_t, 4>> quadNodes;
  std::vector<atlas::gidx_t> quadGidx;

  auto isValidEdge = [&](atlas::gidx_t gidxA, atlas::gidx_t gidxB) -> bool {
    // Convert 1-based globalIdx to 0-based C++ vector index
    const size_t idxA = static_cast<size_t>(gidxA - 1);
    const size_t idxB = static_cast<size_t>(gidxB - 1);

    const double rlonA = lons[idxA] * M_PI / 180.0;
    const double rlatA = lats[idxA] * M_PI / 180.0;
    const double rlonB = lons[idxB] * M_PI / 180.0;
    const double rlatB = lats[idxB] * M_PI / 180.0;

    const double dx = std::cos(rlatA) * std::cos(rlonA) - std::cos(rlatB) * std::cos(rlonB);
    const double dy = std::cos(rlatA) * std::sin(rlonA) - std::cos(rlatB) * std::sin(rlonB);
    const double dz = std::sin(rlatA) - std::sin(rlatB);

    return (dx * dx + dy * dy + dz * dz) <= maxChordSq;
  };

  std::unordered_map<uint64_t, atlas::gidx_t> nodeLookup;  // local node table: (tile, i, j) -> 1-based globalIdx
  nodeLookup.reserve(globalIdx.size());

  for (size_t k = 0; k < globalIdx.size(); ++k) {
    nodeLookup[pack_tile_ij(rank_tiles[k], rank_i[k], rank_j[k])] = globalIdx[k];
  }
  for (size_t n = 0; n < ownedCount_; ++n) {
    int t = rank_tiles[n];
    int i = rank_i[n];
    int j = rank_j[n];

    atlas::gidx_t id_SW = globalIdx[n];

    auto [c_SE, c_NE, c_NW, valid] = get_fv3_quad_corners(t, i, j);
    if (!valid) continue;

    auto it_SE = nodeLookup.find(pack_tile_ij(c_SE.t, c_SE.i, c_SE.j));
    if (it_SE == nodeLookup.end()) continue;

    auto it_NE = nodeLookup.find(pack_tile_ij(c_NE.t, c_NE.i, c_NE.j));
    if (it_NE == nodeLookup.end()) continue;

    auto it_NW = nodeLookup.find(pack_tile_ij(c_NW.t, c_NW.i, c_NW.j));
    if (it_NW == nodeLookup.end()) continue;

    atlas::gidx_t id_SE = it_SE->second;
    atlas::gidx_t id_NE = it_NE->second;
    atlas::gidx_t id_NW = it_NW->second;

   if (!isValidEdge(id_SW, id_SE) ||
        !isValidEdge(id_SE, id_NE) ||
        !isValidEdge(id_NE, id_NW) ||
        !isValidEdge(id_NW, id_SW)) {
      continue;
    }

    quadNodes.push_back({id_SW, id_SE, id_NE, id_NW});
  }

  // Compute globally unique 1-based quad indices across all MPI ranks
  const int localQuads = static_cast<int>(quadNodes.size());
  std::vector<int> quadsPerRank(nproc, 0);
  comm_.allGather(localQuads, quadsPerRank.begin(), quadsPerRank.end());

  atlas::gidx_t quadOffset = 0;
  for (int r = 0; r < my_rank; ++r) {
    quadOffset += quadsPerRank[r];
  }
  quadGidx.resize(localQuads);
  for (int q = 0; q < localQuads; ++q) {
    quadGidx[q] = quadOffset + q + 1;
  }  

/*  compute_quad_nodes(
      ownedCount_,
      N,
      lons, lats,
      rank_tiles, rank_i, rank_j,
      globalIdx,
      quadNodes, quadGidx,
      comm_);
*/

  // only print these if necessary-- later remove
  if (print_debug_info) {
      //Check for duplicate globalIdx in local rank tables
      std::unordered_set<atlas::gidx_t> seen_gidx;
      size_t local_duplicates = 0;
      for (atlas::gidx_t gid : globalIdx) {
        if (!seen_gidx.insert(gid).second) {
          local_duplicates++;
        }
      }
      
      size_t total_duplicates = local_duplicates;
      comm_.allReduceInPlace(total_duplicates, eckit::mpi::Operation::SUM);
      oops::Log::info() << "Total duplicate node registrations across local vectors: " 
                        << total_duplicates << std::endl;
  
      //Validate ghost remoteIdx against owner capacity
      // Share each rank's owned node vector size across all ranks
      std::vector<int> localOwnedCounts(nproc, 0);
      int myOwned = static_cast<int>(ownedCount_);
      comm_.allGather(myOwned, localOwnedCounts.begin(), localOwnedCounts.end());
      
      size_t out_of_bounds_remoteIdx = 0;
      for (size_t k = 0; k < ghosts.size(); ++k) {
        if (ghosts[k] == 1) {
          int owner = partitions[k];
          atlas::idx_t rem = remoteIdx[k];
          // remoteIdx is 1-based, so it must be between 1 and localOwnedCounts[owner]
          if (rem < 1 || rem > localOwnedCounts[owner]) {
            out_of_bounds_remoteIdx++;
          }
        }
      }
      
      size_t total_bad_remotes = out_of_bounds_remoteIdx;
      comm_.allReduceInPlace(total_bad_remotes, eckit::mpi::Operation::SUM);
      oops::Log::info() << "Total out-of-bounds ghost remoteIdx values: " 
                      << total_bad_remotes << std::endl;
  }

  // check lengths are consistent
  const size_t total_nodes = rank_points.size();
  ASSERT(rank_lons.size()   == total_nodes);
  ASSERT(rank_lats.size()   == total_nodes);
  ASSERT(ghosts.size()      == total_nodes);
  ASSERT(globalIdx.size()   == total_nodes);
  ASSERT(remoteIdx.size()   == total_nodes);
  ASSERT(partitions.size()  == total_nodes);
  
  // Check for uninitialized / invalid values
  for (size_t k = 0; k < total_nodes; ++k) {
    ASSERT(globalIdx[k] > 0);
    ASSERT(partitions[k] >= 0 && partitions[k] < comm_.size());
    ASSERT(ghosts[k] == 0 || ghosts[k] == 1);
    if (ghosts[k] == 0) {
      ASSERT(partitions[k] == comm_.rank());
    } else {
      ASSERT(partitions[k] != comm_.rank());
    }
  }

  // 10. Build Atlas Mesh with quad cells 
  eckit::LocalConfiguration meshConf;
  meshConf.set("mpi_comm", comm.name());
  meshConf.set("halo", 0);   // ?? doesn't seem to matter

  const atlas::mesh::MeshBuilder meshBuilder;
  mesh_ = meshBuilder(
      rank_lons, rank_lats, ghosts,
      globalIdx, remoteIdx, /*remote_index_base=*/1, partitions,
      /*triag_nodes_global=*/{}, /*triag_global_index=*/{},
      /*quad_nodes_global=*/quadNodes, /*quad_global_index=*/quadGidx,
      meshConf);
  
  if (print_debug_info) {
      oops::Log::info() << "GeometryLandVector: mesh_nodes_size: " << mesh_.nodes().size()
          << " mesh_cells_size: " << mesh_.cells().size() << std::endl;
  }

  // 11. Global verification that at least one rank holds cells
  size_t local_cells = mesh_.cells().size();
  size_t global_cells = local_cells;
  comm_.allReduceInPlace(global_cells, eckit::mpi::Operation::SUM);
  ASSERT(global_cells > 0);

  // Ensure cells are present 
  ASSERT(mesh_.cells().size() > 0);

  // 12. Instantiate Atlas NodeColumns function space
  functionSpace_ = atlas::functionspace::NodeColumns(mesh_, atlas::util::Config("mpi_comm", comm.name()));  

  lonlat_ = functionSpace_.lonlat();
  ghost_ = functionSpace_.ghost();
  auto glb_idx = functionSpace_.global_index();
  auto r_idx = functionSpace_.remote_index();
  partition_ = functionSpace_.partition();

  auto lonlatView = atlas::array::make_view<double, 2>(lonlat_);
  auto glbidxView = atlas::array::make_view<atlas::gidx_t, 1>(glb_idx);
  auto partView = atlas::array::make_view<int,1>(partition_);
  auto ridxView = atlas::array::make_view<int, 1>(r_idx);
  
  auto ghostView = atlas::array::make_view<int, 1>(ghost_);

  // add global_index field: current name is 'glb_idx
  global_index_ = functionSpace_.createField<atlas::gidx_t>(atlas::option::name("global_index"));
  // Create remote_index field (1-based local index on the remote owning rank)
  remote_index_ = functionSpace_.createField<int>(atlas::option::name("remote_index"));  
  auto global_index_view = atlas::array::make_view<atlas::gidx_t, 1>(global_index_);
  auto remote_index_view  = atlas::array::make_view<int, 1>(remote_index_);
 
  // Explicit 1D latitude and longitude fields
  lat_ = functionSpace_.createField<double>(atlas::option::name("lat") | atlas::option::levels(1));
  lon_ = functionSpace_.createField<double>(atlas::option::name("lon") | atlas::option::levels(1));
  auto latView = atlas::array::make_view<double, 2>(lat_);
  auto lonView = atlas::array::make_view<double, 2>(lon_);

  elevationField_ = functionSpace_.createField<double>(atlas::option::name("height") | atlas::option::levels(1));
  auto elevView   = atlas::array::make_view<double, 2>(elevationField_);
  
  for (size_t i = 0; i < functionSpace_.size(); ++i){
    if (partView(i) == my_rank){  //Populate only owned nodes from the local input arrays
      atlas::gidx_t g_idx = glbidxView(i) - 1;; // 0-based lookup into global elevation array
      elevView(i, 0) = elevations[g_idx];  //TODO: use sorted elevation
      latView(i, 0) = lonlatView(i, 1);
      lonView(i,0) = lonlatView(i, 0);
      global_index_view(i) = glbidxView(i);
      remote_index_view(i) = ridxView(i);
      
    }
  }

  // Create 'owned' field (1 for owned point, 0 for ghost point)
  owned_ = functionSpace_.createField<int>(atlas::option::name("owned") | atlas::option::levels(1));
  auto owned_view = atlas::array::make_view<int, 2>(owned_);
  area_ = functionSpace_.createField<double>(atlas::option::name("area") | atlas::option::levels(1));
  auto area_view = atlas::array::make_view<double, 2>(area_); 
  for (size_t i = 0; i < functionSpace_.size(); ++i) {
      if (partView(i) == my_rank) {
          owned_view(i, 0) = 1;
          area_view(i, 0)  = 1.0;
      } else {
          owned_view(i, 0) = 0;
          area_view(i, 0)  = -1.0;
      }
  }

  // Ensure ghost field matches the manual ghost array
  for (atlas::idx_t j = 0; j < functionSpace_.size(); ++j) {
    int g_expected = (j < ownedCount_) ? 0 : 1;
    ASSERT(ghostView(j) == g_expected);
  }
  
  // Add gmask and mask fields--(expect to use in BUMP, but not sure) 
  atlas::Field gmaskField = functionSpace_.createField<int>(
      atlas::option::name("gmask") | atlas::option::levels(1));
  auto gmaskView = atlas::array::make_view<int, 2>(gmaskField);
  for (atlas::idx_t j = 0; j < functionSpace_.size(); ++j) {
    gmaskView(j, 0) = (j < ownedCount_) ? 1 : 0;
  }
 
  atlas::Field maskField = functionSpace_.createField<int>(
      atlas::option::name("mask") | atlas::option::levels(1));
  auto maskView = atlas::array::make_view<int, 2>(maskField);
  for (atlas::idx_t j = 0; j < functionSpace_.size(); ++j) {
    maskView(j, 0) = (j < ownedCount_) ? 1 : 0;
  }
 
  // do halo exchange
  atlas::functionspace::NodeColumns fsHalo(functionSpace_);

  fsHalo.haloExchange(elevationField_);
  fsHalo.haloExchange(lat_);
  fsHalo.haloExchange(lon_);
  fsHalo.haloExchange(lonlat_);
  //fsHalo.haloExchange(owned_); do not halo exchange owned
  fsHalo.haloExchange(area_);
  fsHalo.haloExchange(global_index_);
  fsHalo.haloExchange(remote_index_);
  fsHalo.haloExchange(gmaskField);
  fsHalo.haloExchange(maskField);

  // add fields to 
  fields_.add(lonlat_);
  fields_.add(global_index_);
  fields_.add(remote_index_);
  fields_.add(ghost_);
  fields_.add(partition_);
  fields_.add(elevationField_);  
  fields_.add(lat_);
  fields_.add(lon_); 
  fields_.add(area_);
  fields_.add(owned_);
  fields_.add(gmaskField);
  fields_.add(maskField);

  functionSpace = functionSpace_;
  geomFields = fields_; 

  // save mesh for inspection
  if (config.has("save_grid_gmsh_to")){

    std::string gFile = config.getString("save_grid_gmsh_to");
    //const atlas::functionspace::NodeColumns fsNC(functionSpace_);
    eckit::LocalConfiguration gmshConf;
    gmshConf.set("coordinates", "xyz");
    gmshConf.set("ghost", true);
    atlas::output::Gmsh gmsh(gFile, gmshConf);
    gmsh.write(mesh_); //fsNC.mesh());

    // only need elevation field for now
    atlas::FieldSet gFields;
    gFields.add(elevationField_);
    gmsh.write(gFields, fsHalo);

    comm_.barrier();
    oops::Log::info() << "GeometryLandVector gmsh saved " << std::endl;
    
    //save to nc too
    std::string grid_resol = config.getString("grid_resol", "C384.mx025");
    std::string ncFile = "landVector_" + grid_resol + "_gFields";
    eckit::LocalConfiguration conf;
    conf.set("filepath", ncFile);
    util::writeFieldSet(comm_, conf, gFields);
    oops::Log::info() << "GeometryLandVector elevation field saved to " << ncFile << std::endl;
  }  

  oops::Log::info() << "GeometryLandVector::GeometryLandVector created" << std::endl;

}

// -----------------------------------------------------------------------------
/*GeometryLandVector::GeometryLandVector(const GeometryLandVector & other)
      :sharedData_(other.sharedData_) { // Copy the shared_ptr (shallow copy
      :comm_(other.comm_),
      functionSpace_(other.functionSpace_),
      fields_(other.fields_),
      lonlat_(other.lonlat_),
      elevationField_(other.elevationField_),
      global_index_(other.global_index_),
      lat_(other.lat_),
      lon_(other.lon_),
      ghost_(other.ghost_),
      partition_(other.partition_),
      remote_index_(other.remote_index_),
      owned_(other.owned_),
      area_(other.area_)  {}*/

// -----------------------------------------------------------------------------
//GeometryLandVector::~GeometryLandVector() {}

// -----------------------------------------------------------------------------
size_t GeometryLandVector::SharedData::globalNodeCount() const {
  /*size_t localSize = functionSpace_.size();
  size_t globalSize = 0;
  // Sum local sizes across all MPI ranks in comm_
  comm_.allReduce(localSize, globalSize, eckit::mpi::Operation::SUM);

  return globalSize;*/
  atlas::functionspace::NodeColumns fsHalo(functionSpace_);
  return fsHalo.nb_nodes_global();
}

// -----------------------------------------------------------------------------
size_t GeometryLandVector::SharedData::localNodeCount() const {
  return functionSpace_.size();
}

// -----------------------------------------------------------------------------
void GeometryLandVector::print(std::ostream & os) const {
  os << "GeometryLandVector: [ Count = " << sharedData_->globalNodeCount() << " points ]";
}

// -----------------------------------------------------------------------------------------------

std::vector<double> GeometryLandVector::verticalCoord(std::string & /*vcUnits*/) const {
  std::stringstream errorMsg;
  errorMsg << "GeometryLandVector::verticalCoord is not implemented" << std::endl;
  throw eckit::NotImplemented(errorMsg.str(), Here());
}

void GeometryLandVector::SharedData::readLatLonElevFromFile(const std::string pathFile, const std::string latName, const std::string lonName, const std::string elevName,
    std::vector<double> & lats, std::vector<double> & lons, std::vector<double> & elevations, 
    std::vector<int>& global_tile, std::vector<int>& global_i, std::vector<int>& global_j)
{
    // NetCDF IDs
    int fileId;

    oops::Log::info() << "Reading file " << pathFile << std::endl;
    nc_rc(nc_open(pathFile.c_str(), NC_NOWRITE, &fileId), "nc_open " + pathFile);

    // Get the variable ID for this field
    int varId;
    nc_rc(nc_inq_varid(fileId, latName.c_str(), &varId), "nc_inq_varid " + latName);

    // Get number of dimensions + their IDs
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    nc_rc(nc_inq_var(fileId, varId,
                     nullptr,   // var name (unused)
                     nullptr,   // type (unused)
                     &ndims,
                     dimids,
                     nullptr),  // attributes (unused)
          "nc_inq_var");

    // Ensure that we are working with a vector (1D variables)
    ASSERT(ndims == 1);
    // Read the variable data
    size_t dimlen;
    nc_rc(nc_inq_dimlen(fileId, dimids[0], &dimlen), "nc_get_dimlen " + latName);
    lats.resize(dimlen);
    lons.resize(dimlen);
    elevations.resize(dimlen);
    global_tile.resize(dimlen);
    global_i.resize(dimlen);
    global_j.resize(dimlen);

    nc_rc(nc_get_var_double(fileId, varId, lats.data()), "nc_get_var_double " + latName);

    nc_rc(nc_inq_varid(fileId, lonName.c_str(), &varId), "nc_inq_varid " + lonName);
    nc_rc(nc_get_var_double(fileId, varId, lons.data()), "nc_get_var_double " + lonName);

    nc_rc(nc_inq_varid(fileId, elevName.c_str(), &varId), "nc_inq_varid " + elevName);
    nc_rc(nc_get_var_double(fileId, varId, elevations.data()), "nc_get_var_double " + elevName);
    
    nc_rc(nc_inq_varid(fileId, "cube_tile", &varId), "nc_inq_varid cube_tile ");  // + elevName);
    nc_rc(nc_get_var_int(fileId, varId, global_tile.data()), "nc_get_var_int cube_tile "); // + elevName);

    nc_rc(nc_inq_varid(fileId, "cube_i", &varId), "nc_inq_varid cube_i" );  
    nc_rc(nc_get_var_int(fileId, varId, global_i.data()), "nc_get_var_int cube_j ");

    nc_rc(nc_inq_varid(fileId, "cube_j", &varId), "nc_inq_varid cube_j ");
    nc_rc(nc_get_var_int(fileId, varId, global_j.data()), "nc_get_var_int cube_j");

    // Close file
    nc_rc(nc_close(fileId), "nc_close "+ pathFile);

    oops::Log::info() << "Done reading file " << pathFile << std::endl;

}

// -----------------------------------------------------------------------------------------------

}  // namespace ijedi

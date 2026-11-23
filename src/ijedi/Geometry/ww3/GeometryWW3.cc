#include "ijedi/Geometry/ww3/GeometryWW3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

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
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"

#include "oops/util/Logger.h"

namespace ijedi
{
namespace
{
// Reads a WW3 gmsh file (v2, ASCII): node lon (0-360), lat and depth (z), and the 0-based
// nodes of the triangles (element type 2; other element types are skipped)
void readMsh(const std::string & path, std::vector<double> & lon, std::vector<double> & lat,
             std::vector<double> & depth, std::vector<int> & tri)
{
  std::ifstream in(path);
  if (!in) throw eckit::CantOpenFile(path, Here());
  std::string line;
  bool format = false, nodes = false;
  while (std::getline(in, line)) {
    if (line.rfind("$MeshFormat", 0) == 0) {
      double version;
      int fileType;
      in >> version >> fileType;
      if (version < 2.0 || version >= 3.0 || fileType != 0)
        throw eckit::BadValue("GeometryWW3: " + path + " is not an ASCII gmsh v2 file", Here());
      format = true;
    } else if (line.rfind("$Nodes", 0) == 0) {
      size_t n;
      in >> n;
      lon.resize(n);
      lat.resize(n);
      depth.resize(n);
      for (size_t i = 0; i < n; ++i) {
        size_t id;
        in >> id >> lon[i] >> lat[i] >> depth[i];
        if (id != i + 1)
          throw eckit::BadValue("GeometryWW3: node ids of " + path + " are not 1..N", Here());
      }
      nodes = true;
    } else if (line.rfind("$Elements", 0) == 0) {
      size_t n;
      in >> n;
      std::getline(in, line);
      for (size_t e = 0; e < n; ++e) {
        std::getline(in, line);
        std::istringstream ss(line);
        int id, type, nTags, tag, n1, n2, n3;
        ss >> id >> type >> nTags;
        if (type != 2) continue;
        for (int t = 0; t < nTags; ++t) ss >> tag;
        ss >> n1 >> n2 >> n3;
        tri.insert(tri.end(), {n1 - 1, n2 - 1, n3 - 1});
      }
    }
  }
  if (!format || !nodes || tri.empty())
    throw eckit::BadValue("GeometryWW3: no format, nodes or triangles in " + path, Here());
}

template <typename T>
void broadcastVector(const eckit::mpi::Comm & comm, std::vector<T> & v)
{
  size_t n = v.size();
  comm.broadcast(n, 0);
  v.resize(n);
  comm.broadcast(v, 0);
}
}  // namespace

// -------------------------------------------------------------------------------------------------
GeometryWW3::GeometryWW3(const eckit::Configuration & conf, const eckit::mpi::Comm & comm,
                         eckit::LocalConfiguration & geomVariables,
                         atlas::FunctionSpace & functionSpace, atlas::FieldSet & geomFields,
                         bool & levelsAreTopDown, int & numLevels)
  : meshFile_(conf.getString("mesh file"))
{
  using atlas::gidx_t;
  using atlas::idx_t;
  const int rank = static_cast<int>(comm.rank());
  const int npes = static_cast<int>(comm.size());

  // 1. Read the msh on rank 0 and broadcast it
  std::vector<double> lon, lat, depth;
  std::vector<int> tri;
  if (rank == 0) readMsh(meshFile_, lon, lat, depth, tri);
  broadcastVector(comm, lon);
  broadcastVector(comm, lat);
  broadcastVector(comm, depth);
  broadcastVector(comm, tri);
  nNodes_ = static_cast<int>(lon.size());
  nTriangles_ = static_cast<int>(tri.size() / 3);

  // 2. equal_regions partition of the nodes. As in GeometryMOM6, every rank computes the same
  //    partition on its own single-task communicator (atlas 32-bit overflow on large grids).
  std::vector<int> partOf(nNodes_);
  {
    std::vector<atlas::PointXY> pts(nNodes_);
    for (int k = 0; k < nNodes_; ++k) pts[k] = atlas::PointXY(lon[k], lat[k]);
    const atlas::UnstructuredGrid ugrid(std::move(pts));
    const std::string selfComm = "ijedi_ww3_partition_self";
    comm.split(rank, selfComm);
    atlas::util::Config partConf;
    partConf.set("type", "equal_regions");
    partConf.set("partitions", npes);
    partConf.set("mpi_comm", selfComm);
    const atlas::grid::Distribution dist = atlas::grid::Partitioner(partConf).partition(ugrid);
    eckit::mpi::deleteComm(selfComm);
    for (int k = 0; k < nNodes_; ++k) partOf[k] = dist.partition(k);
  }

  // 3. Local nodes: owned nodes (in global order) then ghosts. A triangle belongs to the rank
  //    owning its lowest-id node; ghosts are the non-owned nodes of the owned triangles.
  //    Remote index = position on the owning rank, known to every rank from the partition.
  std::vector<int> ownerLocalIdx(nNodes_);
  {
    std::vector<int> count(npes, 0);
    for (int k = 0; k < nNodes_; ++k) ownerLocalIdx[k] = count[partOf[k]]++;
  }
  std::vector<int> local;  // global (0-based) index of each local node
  std::vector<int> localOf(nNodes_, -1);
  for (int k = 0; k < nNodes_; ++k) {
    if (partOf[k] != rank) continue;
    localOf[k] = static_cast<int>(local.size());
    local.push_back(k);
  }
  const int nOwned = static_cast<int>(local.size());

  std::vector<std::array<gidx_t, 3>> triNodes;
  std::vector<gidx_t> triGidx;
  for (int t = 0; t < nTriangles_; ++t) {
    const int * v = &tri[3 * t];
    if (partOf[std::min({v[0], v[1], v[2]})] != rank) continue;
    for (int i = 0; i < 3; ++i) {
      if (localOf[v[i]] < 0) {
        localOf[v[i]] = static_cast<int>(local.size());
        local.push_back(v[i]);
      }
    }
    triNodes.push_back({v[0] + 1, v[1] + 1, v[2] + 1});
    triGidx.push_back(t + 1);
  }
  const int nLocal = static_cast<int>(local.size());

  // 4. Atlas mesh
  std::vector<double> lons(nLocal), lats(nLocal);
  std::vector<int> ghosts(nLocal), partitions(nLocal);
  std::vector<gidx_t> globalIdx(nLocal);
  std::vector<idx_t> remoteIdx(nLocal);
  for (int n = 0; n < nLocal; ++n) {
    const int k = local[n];
    lons[n] = lon[k];
    lats[n] = lat[k];
    ghosts[n] = n < nOwned ? 0 : 1;
    partitions[n] = partOf[k];
    globalIdx[n] = k + 1;
    remoteIdx[n] = ownerLocalIdx[k];
  }
  eckit::LocalConfiguration meshConf;
  meshConf.set("mpi_comm", comm.name());
  atlas::Mesh mesh = atlas::mesh::MeshBuilder()(lons, lats, ghosts, globalIdx, remoteIdx,
                                                 /*remote_index_base=*/0, partitions,
                                                 triNodes, triGidx, {}, {}, meshConf);
  atlas::mesh::actions::build_halo(mesh, 1);
  atlas::functionspace::NodeColumns fs(mesh, atlas::util::Config("mpi_comm", comm.name()));

  // 5. Geometry fields, set on owned nodes and halo-exchanged (NaN check catches any gap)
  atlas::FieldSet fields;
  atlas::Field owned = fs.createField<int>(atlas::option::name("owned")
                                           | atlas::option::levels(1));
  auto vOwned = atlas::array::make_view<int, 2>(owned);
  for (idx_t n = 0; n < owned.shape(0); ++n) vOwned(n, 0) = n < nOwned ? 1 : 0;
  fields.add(owned);
  for (const auto & f : std::vector<std::pair<std::string, const std::vector<double> *>>{
           {"lon", &lon}, {"lat", &lat}, {"depth", &depth}}) {
    atlas::Field field = fs.createField<double>(atlas::option::name(f.first)
                                                | atlas::option::levels(1));
    auto v = atlas::array::make_view<double, 2>(field);
    v.assign(std::numeric_limits<double>::quiet_NaN());
    for (int n = 0; n < nOwned; ++n) v(n, 0) = (*f.second)[local[n]];
    fs.haloExchange(field);
    for (idx_t n = 0; n < v.shape(0); ++n)
      if (std::isnan(v(n, 0)))
        throw eckit::Exception("GeometryWW3: " + f.first + " not set by the halo exchange",
                               Here());
    fields.add(field);
  }

  // 6. Optional gmsh output of the mesh
  if (conf.has("save gmsh to")) {
    eckit::LocalConfiguration gmshConf;
    gmshConf.set("coordinates", conf.getString("gmsh coordinates", "xyz"));
    gmshConf.set("ghost", true);
    atlas::output::Gmsh gmsh(conf.getString("save gmsh to"), gmshConf);
    gmsh.write(mesh);
    atlas::FieldSet toWrite;
    toWrite.add(fields.field("depth"));
    gmsh.write(toWrite, fs);
  }

  // Load-balance statistics for print
  const int nGhost = nLocal - nOwned;
  ownedMin_ = ownedMax_ = nOwned;
  ghostMin_ = ghostMax_ = nGhost;
  comm.allReduceInPlace(ownedMin_, eckit::mpi::min());
  comm.allReduceInPlace(ownedMax_, eckit::mpi::max());
  comm.allReduceInPlace(ghostMin_, eckit::mpi::min());
  comm.allReduceInPlace(ghostMax_, eckit::mpi::max());

  // 7. Outputs
  functionSpace = fs;
  geomFields = fields;
  numLevels = 1;
  levelsAreTopDown = true;
  geomVariables.set("nx", nNodes_);
  geomVariables.set("nLevels", 1);
}

// -------------------------------------------------------------------------------------------------
void GeometryWW3::print(std::ostream & os) const
{
  os << "WW3 unstructured mesh: " << meshFile_ << "\n"
     << "  nodes: " << nNodes_ << ", triangles: " << nTriangles_ << "\n"
     << "  owned nodes per rank: min " << ownedMin_ << " max " << ownedMax_ << "\n"
     << "  ghost nodes per rank: min " << ghostMin_ << " max " << ghostMax_;
}

}  // namespace ijedi

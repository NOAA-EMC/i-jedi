#include <netcdf.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/mpas/GeometryMPAS.h"

namespace
{

  // Thin RAII wrapper around a read-only NetCDF file
  class NcReader
  {
   public:
    explicit NcReader(const std::string &path) : path_(path)
    {
      check(nc_open(path.c_str(), NC_NOWRITE, &ncid_), "nc_open");
    }
    ~NcReader() { nc_close(ncid_); }
    NcReader(const NcReader &) = delete;
    NcReader &operator=(const NcReader &) = delete;

    bool hasDim(const std::string &name) const
    {
      int dimid;
      return nc_inq_dimid(ncid_, name.c_str(), &dimid) == NC_NOERR;
    }

    size_t dimLen(const std::string &name) const
    {
      int dimid;
      size_t len;
      check(nc_inq_dimid(ncid_, name.c_str(), &dimid), "nc_inq_dimid " + name);
      check(nc_inq_dimlen(ncid_, dimid, &len), "nc_inq_dimlen " + name);
      return len;
    }

    bool hasVar(const std::string &name) const
    {
      int varid;
      return nc_inq_varid(ncid_, name.c_str(), &varid) == NC_NOERR;
    }

    // Read rows [row0, row0 + nRows) of a variable whose leading (non-Time) dimension is the
    // row dimension (nCells, nVertices); all trailing dimensions are read whole.
    template <typename T>
    void readRows(const std::string &name, size_t row0, size_t nRows, T *buf) const
    {
      int varid, ndims, unlimdim;
      check(nc_inq_varid(ncid_, name.c_str(), &varid), "nc_inq_varid " + name);
      check(nc_inq_varndims(ncid_, varid, &ndims), "nc_inq_varndims " + name);
      std::vector<int> dimids(ndims);
      check(nc_inq_vardimid(ncid_, varid, dimids.data()), "nc_inq_vardimid " + name);
      check(nc_inq_unlimdim(ncid_, &unlimdim), "nc_inq_unlimdim");

      std::vector<size_t> start(ndims, 0), count(ndims);
      size_t d = 0;
      if (ndims > 0 && dimids[0] == unlimdim) count[d++] = 1;  // first Time level only
      if (d >= static_cast<size_t>(ndims)) {
        throw eckit::BadValue("MPAS variable " + name + " has no row dimension", Here());
      }
      start[d] = row0;
      count[d++] = nRows;
      for (; d < static_cast<size_t>(ndims); ++d) {
        check(nc_inq_dimlen(ncid_, dimids[d], &count[d]), "nc_inq_dimlen " + name);
      }
      getVara(varid, start.data(), count.data(), buf, name);
    }

   private:
    void getVara(int varid, const size_t *start, const size_t *count, int *buf,
                 const std::string &name) const
    {
      check(nc_get_vara_int(ncid_, varid, start, count, buf), "nc_get_vara_int " + name);
    }
    void getVara(int varid, const size_t *start, const size_t *count, double *buf,
                 const std::string &name) const
    {
      check(nc_get_vara_double(ncid_, varid, start, count, buf), "nc_get_vara_double " + name);
    }
    void check(int status, const std::string &what) const
    {
      if (status != NC_NOERR) {
        throw eckit::ReadError("GeometryMPAS: " + what + " failed for " + path_ + ": " +
                               nc_strerror(status), Here());
      }
    }

    std::string path_;
    int ncid_ = -1;
  };

  // -----------------------------------------------------------------------------------------------

  // Read a whole connectivity variable on the root task and broadcast it
  std::vector<int> readAndBroadcast(const NcReader *file, const eckit::mpi::Comm &comm,
                                    const std::string &name, size_t nRows, size_t rowLen)
  {
    std::vector<int> values(nRows * rowLen);
    if (comm.rank() == 0) file->readRows(name, 0, nRows, values.data());
    comm.broadcast(values, 0);
    return values;
  }

  // -----------------------------------------------------------------------------------------------

  // Read rows of a per-cell variable for a list of (0-based) cells. The file is read in
  // contiguous blocks between the first and last requested cell, so each task only reads the
  // span it needs and memory stays bounded however the partition is scattered.
  std::vector<double> readCellRows(const NcReader &file, const std::string &name,
                                   const std::vector<int> &cells, size_t rowLen)
  {
    std::vector<double> out(cells.size() * rowLen);
    if (cells.empty()) return out;

    std::vector<size_t> order(cells.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&cells](size_t a, size_t b) { return cells[a] < cells[b]; });

    const size_t maxBlockValues = 4 * 1024 * 1024;
    const size_t blockRows = std::max<size_t>(1, maxBlockValues / rowLen);
    std::vector<double> block;
    size_t blockStart = 0, blockRowsRead = 0;
    for (const size_t i : order) {
      const size_t cell = cells[i];
      if (blockRowsRead == 0 || cell >= blockStart + blockRowsRead) {
        blockStart = cell;
        blockRowsRead = std::min(blockRows, static_cast<size_t>(cells[order.back()]) + 1 - cell);
        block.resize(blockRowsRead * rowLen);
        file.readRows(name, blockStart, blockRowsRead, block.data());
      }
      std::copy_n(block.begin() + (cell - blockStart) * rowLen, rowLen,
                  out.begin() + i * rowLen);
    }
    return out;
  }

  // -----------------------------------------------------------------------------------------------

  // Owning task of each cell from an MPAS (METIS) graph partition file: one line per cell
  std::vector<int> readPartition(const std::string &filename, const eckit::mpi::Comm &comm,
                                 size_t nCells)
  {
    std::vector<int> owner(nCells);
    if (comm.rank() == 0) {
      std::ifstream in(filename);
      if (!in) {
        throw eckit::CantOpenFile("GeometryMPAS: cannot open partition file " + filename, Here());
      }
      for (size_t c = 0; c < nCells; ++c) {
        if (!(in >> owner[c]) || owner[c] < 0 || owner[c] >= static_cast<int>(comm.size())) {
          throw eckit::BadValue("GeometryMPAS: partition file " + filename + " is not a valid " +
                                std::to_string(comm.size()) + "-way partition of " +
                                std::to_string(nCells) + " cells", Here());
        }
      }
    }
    comm.broadcast(owner, 0);
    return owner;
  }

}  // namespace

namespace ijedi
{

  // -----------------------------------------------------------------------------------------------

  GeometryMPAS::GeometryMPAS(const eckit::Configuration &geomConfig,
                             const eckit::mpi::Comm &comm,
                             eckit::LocalConfiguration &geomVariables,
                             atlas::FunctionSpace &functionSpace,
                             atlas::FieldSet &fieldSet,
                             bool &levelsAreTopDown, int &numberLevels)
  {
    oops::Log::trace() << "GeometryMPAS constructor starting" << std::endl;
    util::Timer timer("ijedi::GeometryMPAS", "GeometryMPAS");

    const int rank = static_cast<int>(comm.rank());
    const std::string gridFile = geomConfig.getString("grid file");
    const int nHalos = geomConfig.getInt("halo layers", 2);  // MPAS config_num_halos
    const std::string vertCoordType = geomConfig.getString("vertical coordinate", "modellevel");
    if (nHalos < 1) {
      throw eckit::BadValue("GeometryMPAS: 'halo layers' must be at least 1", Here());
    }
    if (vertCoordType != "modellevel" && vertCoordType != "height" &&
        vertCoordType != "avgheight") {
      throw eckit::BadValue("GeometryMPAS: unsupported 'vertical coordinate' " + vertCoordType +
                            " (use modellevel, height or avgheight)", Here());
    }

    // Dimensions
    // ----------
    const NcReader file(gridFile);
    const size_t nCells = file.dimLen("nCells");
    const size_t nVertices = file.dimLen("nVertices");
    const size_t maxEdges = file.dimLen("maxEdges");
    const size_t vertexDegree = file.dimLen("vertexDegree");
    if (file.hasDim("nVertLevels")) {
      numberLevels = static_cast<int>(file.dimLen("nVertLevels"));
    } else {
      numberLevels = geomConfig.getInt("number of levels");
    }

    // MPAS stores fields bottom-up (level 1 at the surface)
    levelsAreTopDown = false;

    // Global connectivity, read once and broadcast (freed at the end of the constructor)
    // ----------------------------------------------------------------------------------
    const NcReader *rootFile = rank == 0 ? &file : nullptr;
    const std::vector<int> nEdgesOnCell = readAndBroadcast(rootFile, comm, "nEdgesOnCell",
                                                           nCells, 1);
    const std::vector<int> cellsOnCell = readAndBroadcast(rootFile, comm, "cellsOnCell",
                                                          nCells, maxEdges);
    const std::vector<int> verticesOnCell = readAndBroadcast(rootFile, comm, "verticesOnCell",
                                                             nCells, maxEdges);
    const std::vector<int> cellsOnVertex = readAndBroadcast(rootFile, comm, "cellsOnVertex",
                                                            nVertices, vertexDegree);
    // Regional meshes: the outermost vertices (bdyMaskVertex == 7) form no triangle
    int hasBdyMask = rank == 0 && file.hasVar("bdyMaskVertex") ? 1 : 0;
    comm.broadcast(hasBdyMask, 0);
    const std::vector<int> bdyMaskVertex =
        hasBdyMask ? readAndBroadcast(rootFile, comm, "bdyMaskVertex", nVertices, 1)
                   : std::vector<int>(nVertices, 0);

    // Owner of each cell
    std::vector<int> owner(nCells, 0);
    if (comm.size() > 1) {
      const std::string prefix = geomConfig.getString("partition file prefix");
      owner = readPartition(prefix + std::to_string(comm.size()), comm, nCells);
    }

    // Local cells: owned cells by increasing global index, then each halo layer sorted by
    // global index (as mpas_block_creator_build_cell_halos)
    // ---------------------------------------------------------------------------------
    std::vector<int> localCells;
    std::vector<int> localIndex(nCells, -1);
    for (size_t c = 0; c < nCells; ++c) {
      if (owner[c] == rank) {
        localIndex[c] = static_cast<int>(localCells.size());
        localCells.push_back(static_cast<int>(c));
      }
    }
    const size_t nCellsSolve = localCells.size();
    size_t layerStart = 0;
    for (int layer = 0; layer < nHalos; ++layer) {
      const size_t layerEnd = localCells.size();
      std::vector<int> newCells;
      for (size_t i = layerStart; i < layerEnd; ++i) {
        const size_t c = localCells[i];
        for (int j = 0; j < nEdgesOnCell[c]; ++j) {
          const int n = cellsOnCell[c * maxEdges + j] - 1;
          if (n >= 0 && n < static_cast<int>(nCells) && localIndex[n] < 0) {
            localIndex[n] = 0;  // provisional, so the cell is only added once
            newCells.push_back(n);
          }
        }
      }
      std::sort(newCells.begin(), newCells.end());
      for (const int n : newCells) {
        localIndex[n] = static_cast<int>(localCells.size());
        localCells.push_back(n);
      }
      layerStart = layerEnd;
    }
    const size_t nLocalCells = localCells.size();

    // Index of each local cell on its owning task (owned cells are ordered by global index)
    std::vector<atlas::idx_t> remoteIndex(nLocalCells);
    {
      std::vector<atlas::idx_t> count(comm.size(), 0);
      for (size_t c = 0; c < nCells; ++c) {
        const atlas::idx_t idx = count[owner[c]]++;
        if (localIndex[c] >= 0) remoteIndex[localIndex[c]] = idx;
      }
    }

    // Triangles of the dual mesh, one per owned vertex. Owned vertices are listed as in
    // mpas_block_decomp_all_edges_in_block / _partitioned_edge_list: first encounter over the
    // owned cells' verticesOnCell, owned when the first valid cell on the vertex is owned.
    // --------------------------------------------------------------------------------------
    using atlas::gidx_t;
    std::vector<std::array<gidx_t, 3>> triNodes;
    std::vector<gidx_t> triGlobalIndices;
    if (vertexDegree != 3) {
      throw eckit::BadValue("GeometryMPAS: expected vertexDegree 3 (triangular dual mesh)",
                            Here());
    }
    {
      std::vector<char> seen(nVertices, 0);
      for (size_t i = 0; i < nCellsSolve; ++i) {
        const size_t c = localCells[i];
        for (int j = 0; j < nEdgesOnCell[c]; ++j) {
          const int v = verticesOnCell[c * maxEdges + j] - 1;
          if (v < 0 || seen[v]) continue;
          seen[v] = 1;
          const int *cov = &cellsOnVertex[v * vertexDegree];
          const int *firstValid = std::find_if(cov, cov + vertexDegree,
                                               [](int cell) { return cell > 0; });
          if (owner[*firstValid - 1] != rank || bdyMaskVertex[v] == 7) continue;
          if (std::any_of(cov, cov + vertexDegree, [](int cell) { return cell <= 0; })) continue;
          triNodes.push_back({cov[0], cov[1], cov[2]});  // global cell index = MPAS cell ID
          triGlobalIndices.push_back(v + 1);              // global vertex ID
        }
      }
    }

    // Per-cell coordinates and area for the local cells only
    // ------------------------------------------------------
    const std::vector<double> latCell = readCellRows(file, "latCell", localCells, 1);
    const std::vector<double> lonCell = readCellRows(file, "lonCell", localCells, 1);
    const std::vector<double> areaCell = readCellRows(file, "areaCell", localCells, 1);

    const double rad2deg = 180.0 / M_PI;
    std::vector<double> lons(nLocalCells), lats(nLocalCells);
    std::vector<int> ghosts(nLocalCells), partitions(nLocalCells);
    std::vector<gidx_t> globalIndices(nLocalCells);
    for (size_t i = 0; i < nLocalCells; ++i) {
      lons[i] = lonCell[i] * rad2deg;
      // Single-precision MPAS files put the pole cells slightly beyond +-90 degrees
      lats[i] = std::clamp(latCell[i] * rad2deg, -90.0, 90.0);
      ghosts[i] = i < nCellsSolve ? 0 : 1;
      partitions[i] = owner[localCells[i]];
      globalIndices[i] = localCells[i] + 1;
    }

    // Atlas mesh and function space
    // -----------------------------
    eckit::LocalConfiguration atlasConfig{};
    atlasConfig.set("mpi_comm", comm.name());
    const atlas::mesh::MeshBuilder meshBuilder{};
    atlas::Mesh mesh = meshBuilder(lons, lats, ghosts, globalIndices, remoteIndex, 0, partitions,
                                   triNodes, triGlobalIndices,
                                   std::vector<std::array<gidx_t, 4>>{}, std::vector<gidx_t>{},
                                   atlasConfig);
    atlas::mesh::actions::build_halo(mesh, 1);
    functionSpace = atlas::functionspace::NodeColumns(mesh, atlasConfig);
    const atlas::idx_t fsSize = functionSpace.size();
    ASSERT(fsSize >= static_cast<atlas::idx_t>(nLocalCells));

    if (geomConfig.has("save gmsh to")) {
      atlas::output::Gmsh gmsh(geomConfig.getString("save gmsh to"),
                               atlas::util::Config("coordinates", "xyz") |
                               atlas::util::Config("ghost", true));
      gmsh.write(mesh);
    }

    // Geometry fields: filled on the MPAS local cells, then halo-exchanged so any extra Atlas
    // halo points are filled too
    // ---------------------------------------------------------------------------------------
    fieldSet = atlas::FieldSet();

    atlas::Field owned = functionSpace.createField<int>(atlas::option::name("owned") |
                                                        atlas::option::levels(1));
    auto ownedView = atlas::array::make_view<int, 2>(owned);
    for (atlas::idx_t i = 0; i < fsSize; ++i) {
      ownedView(i, 0) = i < static_cast<atlas::idx_t>(nCellsSolve) ? 1 : 0;
    }
    fieldSet.add(owned);

    atlas::Field area = functionSpace.createField<double>(atlas::option::name("area") |
                                                          atlas::option::levels(1));
    auto areaView = atlas::array::make_view<double, 2>(area);
    for (size_t i = 0; i < nLocalCells; ++i) areaView(i, 0) = areaCell[i];

    // vert_coord: level index ("modellevel"), layer-mid height above sea level ("height") or
    // its global mean over owned cells ("avgheight"), as mpas-jedi's "bump vunit"
    atlas::Field vertCoord = functionSpace.createField<double>(
        atlas::option::name("vert_coord") | atlas::option::levels(numberLevels));
    auto vcView = atlas::array::make_view<double, 2>(vertCoord);
    if (vertCoordType == "modellevel") {
      for (size_t i = 0; i < nLocalCells; ++i) {
        for (int k = 0; k < numberLevels; ++k) vcView(i, k) = k + 1;
      }
    } else {
      const size_t nLevP1 = numberLevels + 1;
      const std::vector<double> zgrid = readCellRows(file, "zgrid", localCells, nLevP1);
      std::vector<double> avgHeight(numberLevels, 0.0);
      for (size_t i = 0; i < nLocalCells; ++i) {
        for (int k = 0; k < numberLevels; ++k) {
          vcView(i, k) = 0.5 * (zgrid[i * nLevP1 + k] + zgrid[i * nLevP1 + k + 1]);
          if (i < nCellsSolve) avgHeight[k] += vcView(i, k);
        }
      }
      if (vertCoordType == "avgheight") {
        comm.allReduceInPlace(avgHeight.begin(), avgHeight.end(), eckit::mpi::sum());
        for (size_t i = 0; i < nLocalCells; ++i) {
          for (int k = 0; k < numberLevels; ++k) vcView(i, k) = avgHeight[k] / nCells;
        }
      }
    }

    atlas::FieldSet haloFields;
    haloFields.add(area);
    haloFields.add(vertCoord);
    functionSpace.haloExchange(haloFields);
    fieldSet.add(area);
    fieldSet.add(vertCoord);

    // Print message
    // -------------
    std::ostringstream msg;
    msg << " MPAS Voronoi mesh from " << gridFile << "\n"
        << " Number of cells: " << nCells << ", vertices: " << nVertices << "\n"
        << " Number of (full) model levels: " << numberLevels << "\n"
        << " Regional mesh: " << (hasBdyMask ? "yes" : "no") << "\n"
        << " Number of tasks: " << comm.size() << ", halo layers: " << nHalos << "\n"
        << " Vertical coordinate: " << vertCoordType;
    printMessage_ = msg.str();

    oops::Log::trace() << "GeometryMPAS constructor done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  void GeometryMPAS::print(std::ostream &os) const
  {
    os << printMessage_ << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  std::vector<double> GeometryMPAS::verticalCoord(std::string &vcUnits) const
  {
    // Not implemented, abort --- IGNORE ---
    std::stringstream errorMsg;
    errorMsg << "GeometryMPAS::verticalCoord is not implemented" << std::endl;
    ABORT(errorMsg.str());
    return {};
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi

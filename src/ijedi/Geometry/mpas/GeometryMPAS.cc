#include <string>
#include <vector>
#include <algorithm>
#include <array>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/option.h"

#include "oops/base/GeometryData.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ijedi/Geometry/mpas/GeometryMPAS.h"
#include "ijedi/Geometry/mpas/GeometryMPAS.interface.h"

namespace {

  // Build an Atlas NodeColumns function space from the MPAS grid topology.
  // Queries node coordinates, ghost/partition metadata, and triangular
  // connectivity from the Fortran geometry object, then assembles the Atlas
  // Mesh and applies a 1-deep halo.
  atlas::functionspace::NodeColumns buildNodeColumns(void * fortranGeom,
                                                     const eckit::mpi::Comm & comm) {
    int num_nodes;
    int num_tri_elements;
    ijedi::ijedi_mpas_geom_get_num_nodes_and_elements_f90(fortranGeom, num_nodes, num_tri_elements);

    std::vector<double> lons(num_nodes);
    std::vector<double> lats(num_nodes);
    std::vector<int> ghosts(num_nodes);
    std::vector<int> global_indices(num_nodes);
    std::vector<int> remote_indices(num_nodes);
    std::vector<int> partitions(num_nodes);

    const int num_tri_nodes = 3 * num_tri_elements;
    std::vector<int> raw_tri_boundary_nodes(num_tri_nodes);
    ijedi::ijedi_mpas_geom_get_coords_and_connectivities_f90(
        fortranGeom, num_nodes, lons.data(), lats.data(), ghosts.data(), global_indices.data(),
        remote_indices.data(), partitions.data(), num_tri_nodes, raw_tri_boundary_nodes.data());

    // Per-PE global tri numbering offset (1-based global element IDs)
    std::vector<int> num_elements_per_rank(comm.size());
    comm.allGather(num_tri_elements, num_elements_per_rank.begin(), num_elements_per_rank.end());
    int global_element_index = 1;
    for (size_t i = 0; i < comm.rank(); ++i) {
      global_element_index += num_elements_per_rank[i];
    }

    using atlas::gidx_t;
    using atlas::idx_t;

    std::vector<gidx_t> atlas_global_indices(num_nodes);
    std::transform(global_indices.begin(), global_indices.end(), atlas_global_indices.begin(),
                   [](const int index) { return atlas::gidx_t{index}; });

    const atlas::idx_t remote_index_base = 1;
    std::vector<idx_t> atlas_remote_indices(num_nodes);
    std::transform(remote_indices.begin(), remote_indices.end(), atlas_remote_indices.begin(),
                   [](const int index) { return atlas::idx_t{index}; });

    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes(num_tri_elements);
    std::vector<gidx_t> tri_global_indices(num_tri_elements);
    for (size_t tri = 0; tri < static_cast<size_t>(num_tri_elements); ++tri) {
      for (size_t i = 0; i < 3; ++i) {
        tri_boundary_nodes[tri][i] = raw_tri_boundary_nodes[3 * tri + i];
      }
      tri_global_indices[tri] = global_element_index++;
    }
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes{};
    std::vector<gidx_t> quad_global_indices{};

    eckit::LocalConfiguration meshConfig{};
    meshConfig.set("mpi_comm", comm.name());
    const atlas::mesh::MeshBuilder mesh_builder{};
    atlas::Mesh mesh = mesh_builder(
        lons, lats, ghosts, atlas_global_indices, atlas_remote_indices, remote_index_base,
        partitions, tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
        quad_global_indices, meshConfig);
    atlas::mesh::actions::build_halo(mesh, 1);
    return atlas::functionspace::NodeColumns(mesh, meshConfig);
  }

  // Fill geometry metadata fields into an Atlas FieldSet using only C++
  // and the geom's scalar metadata + areaCell array exposed via Fortran getters.
  void fillGeometryFields(void * fortranGeom, const atlas::functionspace::NodeColumns & fs,
                          const int nVertLevels, atlas::FieldSet & fields) {
    int nCells, nCellsSolve;
    ijedi::ijedi_mpas_geom_get_local_cell_counts_f90(fortranGeom, nCells, nCellsSolve);

    // All geometry fields (owned, area, vert_coord) loop over fs.size().
    // Assert that Atlas and MPAS agree on the total cell count.
    ASSERT(static_cast<atlas::idx_t>(nCells) == fs.size());

    std::vector<double> areaCell(nCells);
    ijedi::ijedi_mpas_geom_get_area_f90(fortranGeom, nCells, areaCell.data());

    // owned: 1 for solve cells, 0 for halo
    auto owned     = fs.createField<int>(atlas::option::name("owned") | atlas::option::levels(1));
    auto ownedView = atlas::array::make_view<int, 2>(owned);
    for (atlas::idx_t i = 0; i < fs.size(); ++i) {
      ownedView(i, 0) = (i < nCellsSolve) ? 1 : 0;
    }
    fields.add(owned);

    // area
    auto area     = fs.createField<double>(atlas::option::name("area") | atlas::option::levels(1));
    auto areaView = atlas::array::make_view<double, 2>(area);
    for (atlas::idx_t i = 0; i < fs.size(); ++i) {
      areaView(i, 0) = areaCell[i];
    }
    fields.add(area);

    // vert_coord: top-down level index (1 at top, nVertLevels at surface)
    auto vertCoord = fs.createField<double>(atlas::option::name("vert_coord") |
                                            atlas::option::levels(nVertLevels));
    auto vcView    = atlas::array::make_view<double, 2>(vertCoord);
    for (atlas::idx_t i = 0; i < fs.size(); ++i) {
      for (int jz = 0; jz < nVertLevels; ++jz) {
        vcView(i, jz) = static_cast<double>(nVertLevels - jz);
      }
    }
    fields.add(vertCoord);
  }

  std::unordered_map<std::string, size_t> createLevelsPerVariable(const int nVertLevels) {
    return std::unordered_map<std::string, size_t>{
        // JEDI 3D fields on level midpoints
        {"air_horizontal_streamfunction", nVertLevels},
        {"air_horizontal_velocity_potential", nVertLevels},
        {"air_potential_temperature", nVertLevels},
        {"air_pressure", nVertLevels},
        {"air_temperature", nVertLevels},
        {"cloud_liquid_ice", nVertLevels},
        {"cloud_liquid_water", nVertLevels},
        {"dry_air_density", nVertLevels},
        {"eastward_wind", nVertLevels},
        {"effective_radius_of_cloud_ice_particle", nVertLevels},
        {"effective_radius_of_cloud_liquid_water_particle", nVertLevels},
        {"effective_radius_of_graupel_particle", nVertLevels},
        {"effective_radius_of_rain_particle", nVertLevels},
        {"effective_radius_of_snow_particle", nVertLevels},
        {"geopotential_height", nVertLevels},
        {"graupel", nVertLevels},
        {"height_above_mean_sea_level", nVertLevels},
        {"northward_wind", nVertLevels},
        {"mass_content_of_cloud_liquid_water_in_atmosphere_layer", nVertLevels},
        {"mass_content_of_cloud_ice_in_atmosphere_layer", nVertLevels},
        {"mass_content_of_rain_in_atmosphere_layer", nVertLevels},
        {"mass_content_of_snow_in_atmosphere_layer", nVertLevels},
        {"mass_content_of_graupel_in_atmosphere_layer", nVertLevels},
        {"mole_fraction_of_ozone_in_air", nVertLevels},
        {"rain_water", nVertLevels},
        {"relative_humidity", nVertLevels},
        {"snow_water", nVertLevels},
        {"virtual_temperature", nVertLevels},
        {"water_vapor_mixing_ratio_wrt_moist_air", nVertLevels},
        {"water_vapor_mixing_ratio_wrt_dry_air", nVertLevels},
        // JEDI 3D fields on level interfaces
        {"air_pressure_levels", nVertLevels + 1},
        {"geopotential_height_levels", nVertLevels + 1},
        // JEDI 2D fields
        {"air_pressure_at_surface", 1},
        {"average_surface_temperature_within_field_of_view", 1},
        {"eastward_wind_at_surface", 1},
        {"height_above_mean_sea_level_at_surface", 1},
        {"ice_area_fraction", 1},
        {"land_area_fraction", 1},
        {"land_type_index_USGS", 1},
        {"leaf_area_index", 1},
        {"northward_wind_at_surface", 1},
        {"skin_temperature_at_surface_where_ice", 1},
        {"skin_temperature_at_surface_where_land", 1},
        {"skin_temperature_at_surface_where_sea", 1},
        {"skin_temperature_at_surface_where_snow", 1},
        {"soil_temperature", 1},
        {"soil_type", 1},
        {"surface_snow_area_fraction", 1},
        {"surface_snow_thickness", 1},
        {"tropopause_pressure", 1},
        {"vegetation_area_fraction", 1},
        {"vegetation_type_index", 1},
        {"volume_fraction_of_condensed_water_in_soil", 1},
        {"water_area_fraction", 1},
        {"wind_reduction_factor_at_10m", 1}};
  }

}  // namespace

namespace ijedi {

  // -----------------------------------------------------------------------------------------------

  GeometryMPAS::GeometryMPAS(const eckit::Configuration &geomConfig,
                             const eckit::mpi::Comm &comm,
                             eckit::LocalConfiguration &geomVariables,
                             atlas::FunctionSpace &functionSpace,
                             atlas::FieldSet &fieldSet,
                             bool &levelsAreTopDown,
                             int &numberLevels)
  {
  oops::Log::trace() << "GeometryMPAS constructor starting" << std::endl;

  ijedi_mpas_geom_setup_f90(fortranGeom_, geomConfig, &comm);

  // Build Atlas NodeColumns function space from MPAS grid topology
  functionSpace = buildNodeColumns(fortranGeom_, comm);

  // Get nVertLevels from Fortran for vert_coord and levelsPerVariable_
  int nVertLevels;
  ijedi_mpas_geom_get_vertical_resolution_f90(fortranGeom_, nVertLevels);

  // Fill geometry fields (owned, area, vert_coord)
  fieldSet = atlas::FieldSet();
  fillGeometryFields(fortranGeom_, functionSpace, nVertLevels, fieldSet);

  // Set vertical metadata for both the factory outputs and internal storage
  levelsAreTopDown   = false;
  numberLevels       = nVertLevels;
  levelsAreTopDown_  = levelsAreTopDown;
  levelsPerVariable_ = createLevelsPerVariable(numberLevels);

  // Build GeometryData
  geomData_.reset(new oops::GeometryData(functionSpace, fieldSet, levelsAreTopDown, comm));
  comm_ = &comm;  // Store the communicator pointer

  oops::Log::trace() << "ijedi_mpas::GeometryMPAS::GeometryMPAS from config done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  GeometryMPAS::~GeometryMPAS()
  {
    if (fortranGeom_ != nullptr) {
      ijedi_mpas_geom_delete_f90(fortranGeom_);
      fortranGeom_ = nullptr;
    }
  }

  // -----------------------------------------------------------------------------------------------

  void GeometryMPAS::print(std::ostream &os) const
  {
    int nVertLevels;
    ijedi_mpas_geom_get_vertical_resolution_f90(fortranGeom_, nVertLevels);
    int nCellsGlobal;
    ijedi_mpas_geom_get_global_cell_count_f90(fortranGeom_, nCellsGlobal);
    os << "ijedi_mpas::GeometryMPAS, nCellsGlobal = " << nCellsGlobal
       << ", nVertLevels = " << nVertLevels << ", communicator = " << comm_->name();
  }

  // -----------------------------------------------------------------------------------------------

  std::vector<double> GeometryMPAS::verticalCoord(std::string &vcUnits) const
  {
    // Not implemented, abort --- IGNORE ---
    std::stringstream errorMsg;
    errorMsg << "GeometryMPAS::verticalCoord is not implemented" << std::endl;
    ABORT(errorMsg.str());
    return std::vector<double>();  // Never reached, but satisfies compiler
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi

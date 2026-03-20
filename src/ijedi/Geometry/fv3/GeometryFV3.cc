// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/functionspace.h"

#include "ijedi/Geometry/fv3/GeometryFV3.h"

namespace ijedi
{

  GeometryFV3::GeometryFV3(const eckit::Configuration &conf,
                            const eckit::mpi::Comm &comm)
  {
    eckit::mpi::setCommDefault(comm.name().c_str());
    const eckit::LocalConfiguration atlasConfig =
        conf.getSubConfiguration("atlas");

    // Initialize FV3-specific grid dimensions
    // These would typically come from configuration or Fortran backend
    npx_ = conf.getInt("npx", 97);  // Example: C48 resolution has npx=97
    npy_ = conf.getInt("npy", 97);
    npz_ = conf.getInt("npz", 127);  // Vertical levels
    tileNum_ = conf.getInt("tile_num", 1);

    // Set the base class numLevels to match npz
    numLevels_ = npz_;

    // Example: Create Atlas mesh and function space
    // In a real implementation, this would construct an FV3-specific mesh
    // For now, showing the basic pattern:

    // Create a grid (example uses a simple lat-lon grid, FV3 would use
    // cubed-sphere)
    atlas::Grid grid = atlas::Grid(atlasConfig);

    // Generate mesh from grid
    atlas::MeshGenerator meshgen("structured");
    atlas::Mesh mesh = meshgen.generate(grid);

    // Build halos for parallel operations
    // atlas::mesh::actions::build_halo(mesh, 1);

    // Create function space — this is where the fields will be defined
    functionSpace_ = atlas::functionspace::NodeColumns(mesh, atlasConfig);

    // Example: Initialize the fields FieldSet
    // This would contain geometry-related fields like surface height,
    // land mask, etc.
    fields_ = atlas::FieldSet();

    // Example: Add a sample geometry field
    // In a real implementation, these would be read from files or computed
    atlas::Field sampleField = functionSpace_.createField<double>(
        atlas::option::name("surface_altitude") |
        atlas::option::levels(1));
    fields_.add(sampleField);

    // TODO(ijedi): Add more FV3-specific fields as needed:
    // - Land/sea mask
    // - Surface geopotential
    // - Grid cell areas
    // - Coriolis parameter
    // etc.
  }

  void GeometryFV3::print(std::ostream &os) const
  {
    os << "FV3 Geometry:" << std::endl;
    os << "  Resolution (npx x npy): " << npx_ << " x " << npy_
       << std::endl;
    os << "  Number of levels (npz): " << npz_ << std::endl;
    os << "  Tile number: " << tileNum_ << std::endl;
    os << "  Number of fields: " << fields_.size() << std::endl;
  }

  eckit::LocalConfiguration GeometryFV3::gridSpecific() const
  {
    eckit::LocalConfiguration conf;
    conf.set("grid_type", "fv3");
    conf.set("npx", npx_);
    conf.set("npy", npy_);
    conf.set("npz", npz_);
    conf.set("tile_num", tileNum_);
    return conf;
  }

}  // namespace ijedi

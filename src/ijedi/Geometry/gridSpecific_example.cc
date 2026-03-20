// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

// Example of using gridSpecific() to access model-specific parameters
// This demonstrates how State or other classes can safely access grid info

#include <iostream>

#include "ijedi/Geometry/Geometry.h"

namespace ijedi
{
  // Example: How to use gridSpecific in a State class or elsewhere
  void exampleUsage(const Geometry &geom)
  {
    // Get the grid-specific configuration
    eckit::LocalConfiguration gridInfo = geom.gridSpecific();

    // Check what type of grid we're working with
    std::string gridType = gridInfo.getString("grid_type");
    std::cout << "Grid type: " << gridType << std::endl;

    // Safely access parameters based on grid type
    if (gridType == "fv3")
    {
      // For FV3, we can access npx, npy, npz, tile_num
      int npx = gridInfo.getInt("npx");
      int npy = gridInfo.getInt("npy");
      int npz = gridInfo.getInt("npz");
      int tileNum = gridInfo.getInt("tile_num");

      std::cout << "FV3 grid: " << npx << " x " << npy
                << " x " << npz << std::endl;
      std::cout << "Tile number: " << tileNum << std::endl;
    } else if (gridType == "mpas") {
      // For MPAS, we can access nnodes, ncells, nedges
      int nnodes = gridInfo.getInt("nnodes");
      int ncells = gridInfo.getInt("ncells");
      int nedges = gridInfo.getInt("nedges");

      std::cout << "MPAS mesh: " << nnodes << " nodes, "
                << ncells << " cells, " << nedges << " edges"
                << std::endl;
    } else if (gridType == "mom6") {
      // For MOM6, we can access ni, nj, nk
      int ni = gridInfo.getInt("ni");
      int nj = gridInfo.getInt("nj");
      int nk = gridInfo.getInt("nk");

      std::cout << "MOM6 grid: " << ni << " x " << nj
                << " x " << nk << std::endl;
    }

    // Safe access with defaults
    // Returns -999 if not found:
    int someParam = gridInfo.getInt("some_optional_param", -999);

    // Check if a parameter exists before accessing
    if (gridInfo.has("npx"))
    {
      int npx = gridInfo.getInt("npx");
      // Do something FV3-specific...
    }
  }

  // Example: Using in a State constructor
  class StateExample
  {
   public:
    explicit StateExample(const Geometry &geom)
    {
      // Access grid-specific information without knowing the exact type
      auto gridInfo = geom.gridSpecific();

      // Common usage: allocate arrays based on grid dimensions
      std::string gridType = gridInfo.getString("grid_type");
      if (gridType == "fv3")
      {
        int npx = gridInfo.getInt("npx");
        int npy = gridInfo.getInt("npy");
        // Allocate FV3 state arrays...
      }
    }
  };

}  // namespace ijedi

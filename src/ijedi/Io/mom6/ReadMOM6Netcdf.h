#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/field/FieldSet.h"

namespace ijedi {

/// Read a MOM6-format NetCDF file into an Atlas FieldSet.
///
/// Every rank reads the full structured grid from @p filepath, then scatters
/// values into the local (unstructured, ocean-only) Atlas nodes using their
/// global_index.  A halo exchange is performed at the end to fill ghost nodes.
///
/// @param filepath         Input NetCDF path (all ranks read independently).
/// @param x                Distributed FieldSet to populate.
/// @param fileVarNames     Per-field file variable name (empty => skip).
/// @param scalings         Per-field scaling factor (0 => no scaling).
/// @param numLevelsGeom    Number of vertical levels expected by geometry.
void readMOM6Netcdf(const std::string & filepath,
                    atlas::FieldSet & x,
                    const std::vector<std::string> & fileVarNames,
                    const std::vector<double> & scalings,
                    int numLevelsGeom);

}  // namespace ijedi

#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"

#include "eckit/mpi/Comm.h"

namespace ijedi {

/// Read a MOM6-format NetCDF file into an Atlas FieldSet.
///
/// Rank 0 reads the full structured grid from @p filepath, maps it into
/// Atlas global fields using gathered global indices, then scatters to all
/// ranks via Atlas NodeColumns::scatter().  A halo exchange is performed at
/// the end to fill ghost nodes.
///
/// @param filepath         Input NetCDF path (rank 0 reads).
/// @param x                Distributed FieldSet to populate.
/// @param fileVarNames     Per-field file variable name (empty => skip).
/// @param scalings         Per-field scaling factor (0 => no scaling).
/// @param numLevelsGeom    Number of vertical levels expected by geometry.
/// @param comm             MPI communicator.
void readMOM6Netcdf(const std::string & filepath,
                    atlas::FieldSet & x,
                    const std::vector<std::string> & fileVarNames,
                    const std::vector<double> & scalings,
                    int numLevelsGeom,
                    const eckit::mpi::Comm & comm);

}  // namespace ijedi

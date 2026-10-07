#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"

#include "eckit/mpi/Comm.h"

namespace ijedi {

/// Update the variables of an existing CICE restart from an Atlas FieldSet.
///
/// Each field is gathered to rank 0, which overwrites the file variable of the
/// same name at the grid points of the gathered (mesh) nodes; the other grid
/// points and the variables not in the FieldSet keep their values. A field with
/// several levels is written to a variable whose leading dimensions (e.g. ncat)
/// hold that many levels, a single-level field to a (nj, ni) variable.
///
/// @param filepath      Existing NetCDF file, updated on rank 0.
/// @param x             Distributed FieldSet (all ranks participate in gather).
/// @param fileVarNames  Per-field file variable name (empty => skip).
/// @param ni, nj        Global grid dimensions (global index = j * ni + i + 1).
/// @param comm          MPI communicator.
void updateCICENetcdf(const std::string & filepath,
                      const atlas::FieldSet & x,
                      const std::vector<std::string> & fileVarNames,
                      const int ni, const int nj,
                      const eckit::mpi::Comm & comm);

}  // namespace ijedi

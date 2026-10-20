#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"

#include "eckit/mpi/Comm.h"

namespace ijedi {

class GeometryMOM6;

/// Write an Atlas FieldSet to a MOM6-format NetCDF file.
///
/// Each field in @p x is gathered to rank 0 using Atlas NodeColumns::gather(),
/// unpacked into a structured (Time, Layer/lath, lonh) buffer, and written to
/// the file at @p filepath.
///
/// @param filepath      Output NetCDF path (created on rank 0).
/// @param x             Distributed FieldSet (all ranks participate in gather).
/// @param fileVarNames  Per-field file variable name (empty ⇒ skip).
/// @param ni            Global i-dimension (lonh).
/// @param nj            Global j-dimension (lath).
/// @param nz            Number of vertical levels (Layer).
/// @param comm          MPI communicator.
void writeMOM6Netcdf(const std::string & filepath,
                     const atlas::FieldSet & x,
                     const std::vector<std::string> & fileVarNames,
                     int ni, int nj, int nz,
                     const eckit::mpi::Comm & comm);

/// Distributed variant: each task gathers its MOM6 block from the JEDI points
/// and all tasks write their blocks collectively into one NetCDF-4 file
/// (parallel I/O).  No task holds a global field.
void writeMOM6NetcdfBlocks(const std::string & filepath,
                           const atlas::FieldSet & x,
                           const std::vector<std::string> & fileVarNames,
                           int nz,
                           const GeometryMOM6 & geom);

}  // namespace ijedi

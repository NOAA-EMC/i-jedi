#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"

#include "eckit/mpi/Comm.h"

namespace ijedi {

void writeMPASNetcdf(const std::string & filepath,
                     const atlas::FieldSet & x,
                     const std::vector<std::string> & fileVarNames,
                     int nCellsGlobal,
                     int nVertLevels,
                     const eckit::mpi::Comm & comm);

}  // namespace ijedi

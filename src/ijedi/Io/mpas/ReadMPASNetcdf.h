#pragma once

#include <string>
#include <vector>

#include "atlas/field/FieldSet.h"

#include "eckit/mpi/Comm.h"

namespace ijedi {

void readMPASNetcdf(const std::string & filepath,
                    atlas::FieldSet & x,
                    const std::vector<std::string> & fileVarNames,
                    const std::vector<double> & scalings,
                    int nCellsGlobal,
                    const eckit::mpi::Comm & comm);

}  // namespace ijedi

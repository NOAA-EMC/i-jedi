/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Geometry/atlas/GeometryAtlas.h"

#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

namespace ijedi {

// -----------------------------------------------------------------------------------------------

GeometryAtlas::GeometryAtlas(const eckit::Configuration &conf,
                             const eckit::mpi::Comm &comm,
                             eckit::LocalConfiguration & /*geomVariables*/,
                             atlas::FunctionSpace &functionSpace,
                             atlas::FieldSet &geomFields,
                             bool &levelsAreTopDown, int &numLevels)
  : comm_(comm) {
  oops::Log::trace() << "GeometryAtlas constructor starting" << std::endl;

  levelsAreTopDown = conf.getBool("levels are top down", true);
  numLevels_ = conf.getInt("nlevels", 1);
  numLevels = numLevels_;

  // Build the function space from the YAML grid description.
  atlas::Mesh mesh;
  util::setupFunctionSpace(comm, conf, grid_, partitioner_, mesh, functionSpace_, geomFields);
  functionSpace = functionSpace_;

  oops::Log::trace() << "GeometryAtlas constructor done" << std::endl;
}

// -----------------------------------------------------------------------------------------------

void GeometryAtlas::print(std::ostream &os) const {
  std::string prefix;
  os << "Atlas geometry grid:" << std::endl;
  if (grid_)
  {
    os << "- name: " << grid_.name() << std::endl;
    os << "- size: " << grid_.size() << std::endl;
  }
  if (partitioner_)
  {
    os << "Partitioner:" << std::endl;
    os << "- type: " << partitioner_.type() << std::endl;
  }
  os << "Function space:" << std::endl;
  os << "- type: " << functionSpace_.type() << std::endl;
  os << "- levels: " << numLevels_ << std::endl;
}

// -----------------------------------------------------------------------------------------------

std::vector<double> GeometryAtlas::verticalCoord(std::string & /*vcUnits*/) const {
  std::stringstream errorMsg;
  errorMsg << "GeometryAtlas::verticalCoord is not implemented" << std::endl;
  throw eckit::NotImplemented(errorMsg.str(), Here());
}

// -----------------------------------------------------------------------------------------------

}  // namespace ijedi

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
#include "atlas/option.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"

namespace ijedi {

// -----------------------------------------------------------------------------------------------

GeometryAtlas::GeometryAtlas(const eckit::Configuration &conf,
                             const eckit::mpi::Comm &comm,
                             eckit::LocalConfiguration &geomVariables,
                             atlas::FunctionSpace &functionSpace,
                             atlas::FieldSet &geomFields,
                             bool &levelsAreTopDown, int &numLevels)
  : comm_(comm) {
  oops::Log::trace() << "GeometryAtlas constructor starting" << std::endl;

  levelsAreTopDown = conf.getBool("levels are top down", true);

  // An explicit reference pressure column defines the vertical, and with it the level count.
  // "nlevels" stays accepted alongside it, but only as a cross-check: silently disagreeing
  // with the column would misallocate every field in the state.
  if (conf.has("vertical levels")) {
    verticalLevels_ = conf.getDoubleVector("vertical levels");
    if (verticalLevels_.empty()) {
      throw eckit::BadValue("GeometryAtlas: 'vertical levels' must not be empty", Here());
    }
    numLevels_ = static_cast<int>(verticalLevels_.size());
    if (conf.has("nlevels") && conf.getInt("nlevels") != numLevels_) {
      std::stringstream errorMsg;
      errorMsg << "GeometryAtlas: 'nlevels' (" << conf.getInt("nlevels") << ") disagrees with the "
               << "length of 'vertical levels' (" << numLevels_ << ")";
      throw eckit::BadValue(errorMsg.str(), Here());
    }
  } else {
    numLevels_ = conf.getInt("nlevels", 1);
  }
  numLevels = numLevels_;

  // Build the function space from the YAML grid description.
  atlas::Mesh mesh;
  util::setupFunctionSpace(comm, conf, grid_, partitioner_, mesh, functionSpace_, geomFields);
  functionSpace = functionSpace_;

  // Publish what the rest of the stack reads out of the model data.
  geomVariables.set("nLevels", numLevels_);

  if (!verticalLevels_.empty()) {
    // Consumed by ijedi::Geometry under "vertical coordinate source: model", which is the
    // default, and from there by vertical localization.
    geomVariables.set("vertical_coordinate_reference_pressure", verticalLevels_);

    // The same column as a horizontally-uniform 3D field, so that "iterator dimension: 3"
    // (point-mode iteration) has a vertical coordinate to work with.
    atlas::Field vertCoord = functionSpace_.createField<double>(
        atlas::option::name("vert_coord") | atlas::option::levels(numLevels_));
    auto vertCoordView = atlas::array::make_view<double, 2>(vertCoord);
    for (atlas::idx_t j = 0; j < vertCoord.shape(0); ++j) {
      for (int k = 0; k < numLevels_; ++k) {
        vertCoordView(j, k) = verticalLevels_[k];
      }
    }
    geomFields.add(vertCoord);
  }

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
  if (!verticalLevels_.empty()) {
    os << "- vertical coordinate: reference pressure [Pa], "
       << verticalLevels_.front() << " to " << verticalLevels_.back() << std::endl;
  } else {
    os << "- vertical coordinate: none (level count only)" << std::endl;
  }
}

// -----------------------------------------------------------------------------------------------

std::vector<double> GeometryAtlas::verticalCoord(std::string &vcUnits) const {
  if (verticalLevels_.empty()) {
    throw eckit::BadValue("GeometryAtlas::verticalCoord: this geometry was built without "
                          "'vertical levels', so it has no vertical coordinate", Here());
  }
  vcUnits = "Pa";
  return verticalLevels_;
}

// -----------------------------------------------------------------------------------------------

}  // namespace ijedi

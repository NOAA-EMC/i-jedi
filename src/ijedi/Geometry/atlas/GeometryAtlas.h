/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "ijedi/Geometry/base/GeometryBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  // Model-independent geometry that builds an atlas FunctionSpace directly from a YAML grid
  // description.
  //
  // The vertical is described in one of two ways:
  //   - "nlevels": a bare level count, with no coordinate attached (the default).
  //   - "vertical levels": an explicit reference pressure column in Pa, one value per level.
  //     This is what a model discretised on fixed pressure levels needs (NeuralGCM and the
  //     other AI models are all defined on the ERA5 pressure levels). Supplying it also fixes
  //     the level count, publishes the column as "vertical_coordinate_reference_pressure" for
  //     ijedi::Geometry's "vertical coordinate source: model" branch, and builds the
  //     horizontally-uniform "vert_coord" field that point-mode iteration requires.
  class GeometryAtlas : public GeometryBase
  {
   public:
    GeometryAtlas(const eckit::Configuration &, const eckit::mpi::Comm &,
                  eckit::LocalConfiguration &, atlas::FunctionSpace &,
                  atlas::FieldSet &, bool &, int &);
    std::vector<double> verticalCoord(std::string &) const override;

   private:
    void print(std::ostream &) const override;
    const eckit::mpi::Comm & comm_;
    int numLevels_ = 1;
    // Reference pressure per level [Pa]. Empty when the vertical is only a level count.
    std::vector<double> verticalLevels_;
    atlas::Grid grid_;
    atlas::grid::Partitioner partitioner_;
    atlas::FunctionSpace functionSpace_;
  };

}  // namespace ijedi

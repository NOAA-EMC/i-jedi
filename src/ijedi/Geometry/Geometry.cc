// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#include <string>
#include <numeric>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Geometry/base/GeometryBase.h"

// ---------------------------------------------------------------------------
namespace ijedi
{
  // ---------------------------------------------------------------------------
  const int HALO_SIZE = 1;
  // ---------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration &conf,
                     const eckit::mpi::Comm &comm)
      : comm_(comm)
  {
    // Trace
    oops::Log::trace() << "Geometry constructor starting" << std::endl;

    // Create the geometry implementation (which will set numLevels_)
    geometryImpl_ = GeometryBase::create(conf, comm_);

    // Construct the fields metadata object using numLevels from the base
    fieldsMeta_.reset(new FieldsMetadata(geometryImpl_->numLevels()));

    // Trace
    oops::Log::trace() << "Geometry constructor starting" << std::endl;
  }
  // ---------------------------------------------------------------------------
  Geometry::Geometry(const Geometry &other)
      : comm_(other.comm_),
        fieldsMeta_(other.fieldsMeta_),
        geometryImpl_(other.geometryImpl_)
  {
  }
  // ---------------------------------------------------------------------------
  Geometry::~Geometry()
  {
  }
  // ---------------------------------------------------------------------------
  void Geometry::print(std::ostream &os) const
  {
    geometryImpl_->print(os);
  }

  // ---------------------------------------------------------------------------

  std::vector<double> Geometry::verticalCoord(std::string &vcUnits) const
  {
    // Not implemented, abort
    std::stringstream errorMsg;
    errorMsg << "Geometry::verticalCoord is not implemented" << std::endl;
    ABORT(errorMsg.str());

    std::vector<double> vc(geometryImpl_->numLevels());
    return vc;
  }

  // ---------------------------------------------------------------------------

  std::vector<size_t> Geometry::variableSizes(
      const oops::Variables &vars) const
  {
    // Array of level heights
    std::vector<size_t> varSizes;
    // Loop through arrays and search metadata map for the levels
    // for (size_t it = 0; it < vars.size(); it++)
    //{
    //  varSizes.push_back(fieldsMeta_->getLevels(vars[it].name()));
    //}
    return varSizes;
  }

  // ---------------------------------------------------------------------------

}  // namespace ijedi

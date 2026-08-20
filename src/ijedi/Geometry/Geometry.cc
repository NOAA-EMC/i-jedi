#include <string>
#include <numeric>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Geometry/base/GeometryBase.h"
#include "mist/utils/Geometry.h"
#include "mist/utils/IterationUnit.h"

// -------------------------------------------------------------------------------------------------
namespace ijedi
{
  // -----------------------------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration &geomConf, const eckit::mpi::Comm &comm)
      : mist::Geometry(comm)
  {
    // Trace
    oops::Log::trace() << "Geometry constructor starting" << std::endl;

    // Get the type
    if (geomConf.has("geometry_type"))
    {
      type_ = geomConf.getString("geometry_type");
    } else {
      // Abort
      std::stringstream errorMsg;
      errorMsg << "Geometry type (geometry_type) not specified in configuration.";
      throw eckit::BadValue(errorMsg.str(), Here());
    }

    // Create the geometry implementation (which will set numLevels_)
    geometryImpl_ = GeometryBase::create(geomConf, comm, modelData_, functionspace_,
                                         fields_, levelsAreTopDown_, numberLevels_);

    // Construct the fields metadata object using numLevels from the base class
    fieldsMeta_ = std::make_shared<FieldsMetadata>(numberLevels_);

    // Set up levels information for each variable using the fields metadata
    levelsPerVariable_ = fieldsMeta_->levelsPerVariable();

    // Expose vertical ordering to downstream components such as Vader recipes.
    modelData_.set("levels_are_top_down", levelsAreTopDown_);

    // Reference pressure column used by vertical localization. Prefer one published by
    // the model geometry (which can follow that model's own mid-layer convention, see
    // GeometryFV3), then mist's generic hybrid-sigma version, then level indices.
    if (modelData_.has("vertical_coordinate_reference_pressure")) {
      verticalCoord_ = modelData_.getDoubleVector("vertical_coordinate_reference_pressure");
    } else if (modelData_.has("sigma_pressure_hybrid_coordinate_a_coefficient") &&
               modelData_.has("sigma_pressure_hybrid_coordinate_b_coefficient")) {
      const std::vector<double> ak =
          modelData_.getDoubleVector("sigma_pressure_hybrid_coordinate_a_coefficient");
      const std::vector<double> bk =
          modelData_.getDoubleVector("sigma_pressure_hybrid_coordinate_b_coefficient");
      verticalCoord_ = std::get<0>(mist::setupReferencePressure(functionspace_, ak, bk));
    } else {
      verticalCoord_.resize(numberLevels_);
      std::iota(verticalCoord_.begin(), verticalCoord_.end(), 0.0);
    }

    // Build GeometryData
    geomData_.reset(new oops::GeometryData(functionspace_, fields_, levelsAreTopDown_, comm));

    // Enable grid-point iteration. Default of 2 (whole columns) is what
    // oops::VerticalLocEV relies on when populating its eigenvectors.
    const int iteratorDimension = geomConf.getInt("iterator dimension", 2);
    if (iteratorDimension != 2 && iteratorDimension != 3) {
      throw eckit::BadValue("ijedi::Geometry: 'iterator dimension' must be 2 or 3", Here());
    }
    initIteratorSupport(iteratorDimension == 2 ? mist::IterationUnit::Column
                                                : mist::IterationUnit::Point);

    // Trace
    oops::Log::trace() << "Geometry constructor finished" << std::endl;
  }
  // -----------------------------------------------------------------------------------------------
  Geometry::~Geometry()
  {
  }
  // -----------------------------------------------------------------------------------------------
  void Geometry::print(std::ostream &os) const
  {
    // Write a general message about the geometry and the implementation provider
    os << std::endl
       << "--------------------------------------------------"
          "--------------------------------------------------";
    os << std::endl
       << "Geometry Information (from type: " + type_ + "):" << std::endl;
    geometryImpl_->print(os);
    os << std::endl
       << "--------------------------------------------------"
          "--------------------------------------------------";
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi

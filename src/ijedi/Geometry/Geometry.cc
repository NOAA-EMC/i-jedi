#include <string>
#include <numeric>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Geometry/base/GeometryBase.h"

// -------------------------------------------------------------------------------------------------
namespace ijedi
{
  // -----------------------------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration &geomConf, const eckit::mpi::Comm &comm)
      : mist::base::Geometry(comm)
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

    // Build GeometryData
    geomData_.reset(new oops::GeometryData(functionspace_, fields_, levelsAreTopDown_, comm));

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

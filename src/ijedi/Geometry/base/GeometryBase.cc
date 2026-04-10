#include <algorithm>
#include <memory>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "ijedi/Geometry/base/GeometryBase.h"
#include "ijedi/Geometry/fv3/GeometryFV3.h"
#include "ijedi/Geometry/mom6/GeometryMOM6.h"

namespace ijedi
{
  // Implemented in GeometryMPASLoader.cc (compiled into libijedi.so).
  // Loads libijedi_mpas.so on demand and delegates construction to it.
  std::shared_ptr<GeometryBase> loadMPASGeometry(const eckit::Configuration &,
                                                  const eckit::mpi::Comm &,
                                                  eckit::LocalConfiguration &,
                                                  atlas::FunctionSpace &,
                                                  atlas::FieldSet &,
                                                  bool &, int &);

  std::shared_ptr<GeometryBase> GeometryBase::create(const eckit::Configuration &geomConf,
                                                     const eckit::mpi::Comm &comm,
                                                     eckit::LocalConfiguration &geomVars,
                                                     atlas::FunctionSpace &functionSpace,
                                                     atlas::FieldSet &fieldSet,
                                                     bool &levelsAreTopDown, int &numLevels)
  {
    // Get the type
    std::string type;
    type = geomConf.getString("geometry_type");

    if (type == "fv3")
    {
      return std::make_shared<GeometryFV3>(geomConf, comm, geomVars, functionSpace, fieldSet,
                                           levelsAreTopDown, numLevels);
    }
    if (type == "mpas")
    {
      return loadMPASGeometry(geomConf, comm, geomVars, functionSpace, fieldSet,
                               levelsAreTopDown, numLevels);
    }
    if (type == "mom6")
    {
      return std::make_shared<GeometryMOM6>(geomConf, comm, geomVars, functionSpace, fieldSet,
                                            levelsAreTopDown, numLevels);
    }

    throw eckit::BadValue("Unsupported geometry type: " + type,
                          Here());
  }

}  // namespace ijedi

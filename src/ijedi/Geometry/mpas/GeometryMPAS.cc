#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ijedi/Geometry/mpas/GeometryMPAS.h"

namespace ijedi
{

  // -----------------------------------------------------------------------------------------------

  GeometryMPAS::GeometryMPAS(const eckit::Configuration &geomConfig,
                             const eckit::mpi::Comm &comm,
                             eckit::Configuration &geomVariables,
                             atlas::FunctionSpace &functionSpace,
                             atlas::FieldSet &fieldSet,
                             int &numberLevels) {}

  // -----------------------------------------------------------------------------------------------

  void GeometryMPAS::print(std::ostream &os) const {}

  // -----------------------------------------------------------------------------------------------

  std::vector<double> GeometryMPAS::verticalCoord(std::string &vcUnits) const
  {
    // Not implemented, abort --- IGNORE ---
    std::stringstream errorMsg;
    errorMsg << "GeometryMPAS::verticalCoord is not implemented" << std::endl;
    ABORT(errorMsg.str());
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi

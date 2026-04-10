#include <cstdlib>

#include "eckit/config/Configuration.h"

#include "ijedi/Geometry/base/GeometryBase.h"
#include "ijedi/Geometry/mpas/GeometryMPAS.h"
#include "ijedi/Geometry/mpas/GeometryMPAS.interface.h"

namespace ijedi {

namespace {

void mpasGeomESMFAtExit() { ijedi_mpas_esmf_shutdown_f90(); }
bool registerESMFAtExit() { std::atexit(mpasGeomESMFAtExit); return true; }

}  // namespace

extern "C" GeometryBase * ijedi_create_geometry_mpas(const eckit::Configuration & geomConf,
                                                       const eckit::mpi::Comm & comm,
                                                       eckit::Configuration & geomVars,
                                                       atlas::FunctionSpace & functionSpace,
                                                       atlas::FieldSet & fieldSet,
                                                       bool & levelsAreTopDown,
                                                       int & numLevels) {
  static bool registered = registerESMFAtExit();
  (void)registered;

  return new GeometryMPAS(geomConf, comm, geomVars, functionSpace, fieldSet,
                          levelsAreTopDown, numLevels);
}

}  // namespace ijedi

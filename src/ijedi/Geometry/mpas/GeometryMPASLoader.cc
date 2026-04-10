#include <dlfcn.h>
#include <sstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "ijedi/Geometry/base/GeometryBase.h"

namespace ijedi
{
  namespace {

    using CreateGeometryMPASFn = GeometryBase * (*)(const eckit::Configuration &,
                                                    const eckit::mpi::Comm &,
                                                    eckit::LocalConfiguration &,
                                                    atlas::FunctionSpace &,
                                                    atlas::FieldSet &,
                                                    bool &, int &);

    CreateGeometryMPASFn loadMPASGeometryFactory(std::string &errMsg)
    {
      static void * mpasPluginHandle = nullptr;
      static CreateGeometryMPASFn createFn = nullptr;
      static bool attemptedLoad = false;

      if (createFn != nullptr) {
        return createFn;
      }

      if (attemptedLoad) {
        errMsg = "MPAS geometry plugin has been requested before but could not be loaded.";
        return nullptr;
      }

      attemptedLoad = true;

      dlerror();
      mpasPluginHandle = dlopen("libijedi_mpas.so", RTLD_NOW | RTLD_LOCAL);
      if (mpasPluginHandle == nullptr) {
        const char * loadErr = dlerror();
        errMsg = loadErr ? loadErr : "unknown error while loading libijedi_mpas.so";
        return nullptr;
      }

      dlerror();
      void * symbol = dlsym(mpasPluginHandle, "ijedi_create_geometry_mpas");
      if (symbol == nullptr) {
        const char * symErr = dlerror();
        errMsg = symErr ? symErr : "symbol ijedi_create_geometry_mpas not found";
        return nullptr;
      }

      createFn = reinterpret_cast<CreateGeometryMPASFn>(symbol);
      return createFn;
    }

  }  // namespace

  std::shared_ptr<GeometryBase> loadMPASGeometry(const eckit::Configuration &geomConf,
                                                  const eckit::mpi::Comm &comm,
                                                  eckit::LocalConfiguration &geomVars,
                                                  atlas::FunctionSpace &functionSpace,
                                                  atlas::FieldSet &fieldSet,
                                                  bool &levelsAreTopDown, int &numLevels)
  {
    std::string dlErr;
    auto createMPAS = loadMPASGeometryFactory(dlErr);
    if (createMPAS == nullptr) {
      std::stringstream errorMsg;
      errorMsg << "MPAS geometry requested, but MPAS plugin could not be loaded: " << dlErr;
      throw eckit::BadValue(errorMsg.str(), Here());
    }

    return std::shared_ptr<GeometryBase>(
        createMPAS(geomConf, comm, geomVars, functionSpace, fieldSet, levelsAreTopDown, numLevels));
  }

}  // namespace ijedi

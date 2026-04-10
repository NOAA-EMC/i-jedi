#include <algorithm>
#include <dlfcn.h>
#include <memory>
#include <string>
#include <sstream>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "ijedi/Geometry/base/GeometryBase.h"
#include "ijedi/Geometry/fv3/GeometryFV3.h"
#include "ijedi/Geometry/mom6/GeometryMOM6.h"

namespace ijedi
{
  namespace {

    using CreateGeometryMPASFn = GeometryBase * (*)(const eckit::Configuration &,
                                                    const eckit::mpi::Comm &,
                                                    eckit::Configuration &,
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

  std::shared_ptr<GeometryBase> GeometryBase::create(const eckit::Configuration &geomConf,
                                                     const eckit::mpi::Comm &comm,
                                                     eckit::Configuration &geomVars,
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
      std::string dlErr;
      CreateGeometryMPASFn createMPAS = loadMPASGeometryFactory(dlErr);
      if (createMPAS == nullptr) {
        std::stringstream errorMsg;
        errorMsg << "MPAS geometry requested, but MPAS plugin could not be loaded: " << dlErr;
        throw eckit::BadValue(errorMsg.str(), Here());
      }

      return std::shared_ptr<GeometryBase>(
          createMPAS(geomConf, comm, geomVars, functionSpace, fieldSet, levelsAreTopDown,
                     numLevels));
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

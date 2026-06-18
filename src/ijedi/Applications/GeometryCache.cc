#include "ijedi/Applications/GeometryCache.h"

#include <ostream>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ijedi/Geometry/Geometry.h"

namespace ijedi {

GeometryCache::GeometryCache(const eckit::mpi::Comm & comm) : Application(comm) {}

int GeometryCache::execute(const eckit::Configuration & fullConfig) const {
  if (!fullConfig.has("geometry")) {
    throw eckit::BadParameter("geometry_cache requires a 'geometry' section in the config", Here());
  }

  const eckit::LocalConfiguration geometryConfig(fullConfig, "geometry");
  const Geometry geometry(geometryConfig, this->getComm());

  oops::Log::info() << "geometry_cache: geometry constructed successfully" << std::endl;
  oops::Log::test() << geometry << std::endl;

  return 0;
}

std::string GeometryCache::appname() const {
  return "ijedi::GeometryCache";
}

}  // namespace ijedi

#pragma once

#include <string>

#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

namespace eckit {
  class Configuration;
}

namespace eckit {
namespace mpi {
  class Comm;
}  // namespace mpi
}  // namespace eckit

namespace ijedi {

class GeometryCache : public oops::Application {
 public:
  explicit GeometryCache(const eckit::mpi::Comm & comm = oops::mpi::world());

  int execute(const eckit::Configuration & fullConfig) const override;

 private:
  std::string appname() const override;
};

}  // namespace ijedi

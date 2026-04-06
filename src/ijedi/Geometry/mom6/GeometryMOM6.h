#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "ijedi/Geometry/base/GeometryBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class GeometryMOM6 : public GeometryBase
  {
   public:
    GeometryMOM6(const eckit::Configuration &, const eckit::mpi::Comm &, eckit::Configuration &,
                 atlas::FunctionSpace &, atlas::FieldSet &, int &);
    void print(std::ostream &) const override;
    std::vector<double> verticalCoord(std::string &) const override;
  };

}  // namespace ijedi

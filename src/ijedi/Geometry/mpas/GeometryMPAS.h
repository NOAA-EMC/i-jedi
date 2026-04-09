#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "atlas/grid.h"

#include "ijedi/Geometry/base/GeometryBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class GeometryMPAS : public GeometryBase
  {
   public:
    GeometryMPAS(const eckit::Configuration &, const eckit::mpi::Comm &,
                 eckit::LocalConfiguration &,
                 atlas::FunctionSpace &, atlas::FieldSet &, bool &, int &);
    void print(std::ostream &) const override;
    std::vector<double> verticalCoord(std::string &) const override;
  };

}  // namespace ijedi

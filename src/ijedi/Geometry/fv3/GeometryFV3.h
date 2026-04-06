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

  class GeometryFV3 : public GeometryBase
  {
   public:
    GeometryFV3(const eckit::Configuration &, const eckit::mpi::Comm &,
                eckit::Configuration &, atlas::FunctionSpace &, atlas::FieldSet &, int &);
    void print(std::ostream &) const override;
    std::vector<double> verticalCoord(std::string &) const override;

   private:
    std::string printMessage_;
  };

}  // namespace ijedi

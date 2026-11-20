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

  // MPAS Voronoi mesh read directly from an MPAS NetCDF file (no MPAS-Model library). The
  // decomposition follows the MPAS framework so the function space matches mpas-jedi.
  class GeometryMPAS : public GeometryBase
  {
   public:
    GeometryMPAS(const eckit::Configuration &, const eckit::mpi::Comm &,
                 eckit::LocalConfiguration &,
                 atlas::FunctionSpace &, atlas::FieldSet &, bool &, int &);
    void print(std::ostream &) const override;
    std::vector<double> verticalCoord(std::string &) const override;

   private:
    std::string printMessage_;
  };

}  // namespace ijedi

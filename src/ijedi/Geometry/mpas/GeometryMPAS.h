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
    GeometryMPAS(const eckit::Configuration &, const eckit::mpi::Comm &, eckit::Configuration &,
                 atlas::FunctionSpace &, atlas::FieldSet &, int &);
    void print(std::ostream &) const override;
    std::vector<double> verticalCoord(std::string &) const override;

    GeometryMPAS(const eckit::Configuration & config, const eckit::mpi::Comm & comm);
    GeometryMPAS(const GeometryMPAS & other);
    ~GeometryMPAS();

    GeometryMPAS & operator=(const GeometryMPAS &) = delete;

    void * fortranGeom() const { return fortranGeom_; }

   private:

    void * fortranGeom_ = nullptr;

  };

}  // namespace ijedi

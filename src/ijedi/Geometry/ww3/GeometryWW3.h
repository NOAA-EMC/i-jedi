#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "ijedi/Geometry/base/GeometryBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  /// WW3 unstructured grid, built from the gmsh (v2, ASCII) file WW3 itself reads.
  ///
  /// The mesh keeps WW3's nodes and triangles. Node global index = msh node id, which is
  /// also the node index (nx) of the WW3 netCDF output. Every node is a wet, computed node:
  /// the geometry has no mask, so fields are interpolated unmasked.
  class GeometryWW3 : public GeometryBase
  {
   public:
    GeometryWW3(const eckit::Configuration &, const eckit::mpi::Comm &,
                eckit::LocalConfiguration &, atlas::FunctionSpace &,
                atlas::FieldSet &, bool &, int &);
    void print(std::ostream &) const override;
    std::vector<double> verticalCoord(std::string &) const override { return {1.0}; }

   private:
    std::string meshFile_;
    int nNodes_ = 0;
    int nTriangles_ = 0;
    int ownedMin_ = 0, ownedMax_ = 0;
    int ghostMin_ = 0, ghostMax_ = 0;
  };

}  // namespace ijedi

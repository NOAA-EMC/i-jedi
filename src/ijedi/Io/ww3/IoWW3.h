#pragma once

#include <ostream>
#include <string>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/IoBase.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  class IoWW3Parameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoWW3Parameters, IoParametersBase)

   public:
    oops::Parameter<std::string> datapath{"datapath", "path to location of files", "./", this};
    oops::RequiredParameter<std::string> filename{"filename", "netCDF file name", this};
    oops::Parameter<std::string> nodeDimension{"node dimension", "dimension of the nodes",
                                               "nx", this};
  };

  // -------------------------------------------------------------------------------------------------

  /// Io for 2-D fields on the WW3 unstructured grid (geometry_type: ww3).
  ///
  /// read: variables on the node dimension (nx in a WW3 restart or field output; other dimensions
  /// of length 1), in msh node order. _FillValue becomes the missing value.
  /// write: a new file with the node dimension, one float variable per field, and lon/lat.
  /// Both use the generic readNetcdf / updateNetcdf (Io/utils/Netcdf.h).
  class IoWW3 : public IoBase, private util::ObjectCounter<IoWW3>
  {
   public:
    static const std::string classname() { return "ijedi::IoWW3"; }

    typedef IoWW3Parameters Parameters_;

    IoWW3(const Geometry &, const Parameters_ &);
    ~IoWW3() = default;
    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry & geom_;
    std::string filepath_;
    std::string nodeDim_;
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

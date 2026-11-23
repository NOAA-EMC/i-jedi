#include "ijedi/Io/ww3/IoWW3.h"

#include <netcdf.h>

#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/utils/Netcdf.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoWW3> makerIoWW3_("ww3");
  // -------------------------------------------------------------------------------------------------
  IoWW3::IoWW3(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()), geom_(geom),
        filepath_(params.datapath.value() + "/" + params.filename.value()),
        nodeDim_(params.nodeDimension.value())
  {}

  // -------------------------------------------------------------------------------------------------
  void IoWW3::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                   const eckit::LocalConfiguration &) const
  {
    util::Timer timer(classname(), "read");
    std::vector<std::string> names;
    for (const auto & field : x) names.push_back(fileionames.getString(field.name(),
                                                                       field.name()));
    readNetcdf(filepath_, x, names, NetcdfLayout{{nodeDim_}, {}}, geom_.getComm());
  }

  // -------------------------------------------------------------------------------------------------
  void IoWW3::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                    const eckit::LocalConfiguration &) const
  {
    util::Timer timer(classname(), "write");
    const eckit::mpi::Comm & comm = geom_.getComm();

    // The fields and the node coordinates
    atlas::FieldSet out;
    std::vector<std::string> names;
    for (const auto & field : x) {
      out.add(field);
      names.push_back(fileionames.getString(field.name(), field.name()));
    }
    for (const std::string name : {"lon", "lat"}) {
      out.add(geom_.fields().field(name));
      names.push_back(name);
    }

    // Rank 0 creates the file: the node dimension and one float variable per field, whose
    // _FillValue is the missing value, then updateNetcdf writes the values. The double missing
    // value is a float value, so missing values are exactly the fill value in the file.
    int status = NC_NOERR;
    if (comm.rank() == 0) {
      const float fill = static_cast<float>(util::missingValue<double>());
      int ncid, dimid, varid;
      status = nc_create(filepath_.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid);
      if (status == NC_NOERR) {
        status = nc_def_dim(ncid, nodeDim_.c_str(), geom_.modelData().getInt("nx"), &dimid);
        for (const auto & name : names) {
          if (status == NC_NOERR) status = nc_def_var(ncid, name.c_str(), NC_FLOAT, 1, &dimid,
                                                      &varid);
          if (status == NC_NOERR) status = nc_put_att_float(ncid, varid, "_FillValue", NC_FLOAT,
                                                            1, &fill);
        }
        const int closeStatus = nc_close(ncid);
        if (status == NC_NOERR) status = closeStatus;
      }
    }
    comm.broadcast(status, 0);
    if (status != NC_NOERR)
      throw eckit::UserError(classname() + ": cannot create " + filepath_ + ": "
                             + nc_strerror(status), Here());
    updateNetcdf(filepath_, out, names, NetcdfLayout{{nodeDim_}, {}}, comm);
  }

  // -------------------------------------------------------------------------------------------------
  void IoWW3::print(std::ostream &os) const
  {
    os << classname() << " Io for 2-D fields on the WW3 grid: " << filepath_;
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

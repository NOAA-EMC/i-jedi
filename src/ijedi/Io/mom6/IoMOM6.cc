#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/mom6/IoMOM6.h"
#include "ijedi/Io/mom6/ReadMOM6Netcdf.h"
#include "ijedi/Io/mom6/WriteMOM6Netcdf.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoMOM6> makerIoMOM6_("mom6");
  // -------------------------------------------------------------------------------------------------
  IoMOM6::IoMOM6(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()),
        geom_(geom),
        filename_(params.filename.value().value_or("MOM.res.nc"))
  {
    util::Timer timer(classname(), "IoMOM6");
    oops::Log::trace() << classname() << " constructor starting" << std::endl;
    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoMOM6::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                    const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read state");
    oops::Log::trace() << classname() << " read state starting" << std::endl;

    // Get geometry info for level count check
    const auto modelData = geom_.modelData();
    const int nz = modelData.getInt("nz");

    // Resolve variable names and scalings for all fields
    std::vector<std::string> fileVarNames;
    std::vector<double> scalings;
    for (const auto & field : x) {
      const std::string jediName = field.name();
      fileVarNames.push_back(fileionames.has(jediName)
                             ? fileionames.getString(jediName) : "");
      scalings.push_back(fileioscaling.has(jediName)
                         ? fileioscaling.getDouble(jediName) : 0.0);
    }

    // Delegate to the NetCDF reader (root reads, scatter, halo exchange)
    readMOM6Netcdf(filename_, x, fileVarNames, scalings, nz,
                   geom_.getComm());

    oops::Log::trace() << classname() << " read state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoMOM6::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                     const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;

    // Get geometry info
    const auto modelData = geom_.modelData();
    const int ni = modelData.getInt("ni");
    const int nj = modelData.getInt("nj");
    const int nz = modelData.getInt("nz");

    // Resolve variable names for all fields
    std::vector<std::string> fileVarNames;
    for (const auto & field : x) {
      const std::string jediName = field.name();
      fileVarNames.push_back(fileionames.has(jediName)
                             ? fileionames.getString(jediName) : "");
    }

    // Delegate to the NetCDF writer (Atlas gather + root writes)
    writeMOM6Netcdf(filename_, x, fileVarNames, ni, nj, nz, geom_.getComm());

    oops::Log::trace() << classname() << " write state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoMOM6::print(std::ostream &os) const
  {
    os << classname() << " Io for MOM6 restarts and histories";
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

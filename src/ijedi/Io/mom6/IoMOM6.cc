#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/mom6/FillMaskedCells.h"
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
        geom_(geom)
  {
    util::Timer timer(classname(), "IoMOM6");
    oops::Log::trace() << classname() << " constructor starting" << std::endl;

    const std::string datapath = params.datapath.value();

    writeFilepath_ = datapath + "/" + params.ocn_file.value();
    readFilepaths_.push_back(writeFilepath_);

    const auto addFile = [&](const boost::optional<std::string> & opt) {
      if (opt) readFilepaths_.push_back(datapath + "/" + *opt);
    };
    addFile(params.ice_file.value());
    addFile(params.fix_file.value());

    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoMOM6::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                    const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read state");
    oops::Log::trace() << classname() << " read state starting" << std::endl;

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
    readMOM6Netcdf(readFilepaths_, x, fileVarNames, scalings, geom_.getComm());

    // Apply boundary conditions to masked (land) cells: nearest-neighbour
    // extrapolation for tracers, zero for non-tracers.
    applyBoundaryConditions(x, geom_.fields().field("mask3d"),
                          geom_.getFieldMetadata());

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
    writeMOM6Netcdf(writeFilepath_, x, fileVarNames, ni, nj, nz, geom_.getComm());

    oops::Log::trace() << classname() << " write state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoMOM6::print(std::ostream &os) const
  {
    os << classname() << " Io for MOM6 ocean, sea ice, and fix files";
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

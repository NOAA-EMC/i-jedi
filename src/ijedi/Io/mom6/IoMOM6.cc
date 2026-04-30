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
        geom_(geom)
  {
    util::Timer timer(classname(), "IoMOM6");
    oops::Log::trace() << classname() << " constructor starting" << std::endl;

    const std::string datapath = params.datapath.value();

    // Build the ordered list of files to read from named parameters
    const auto addFile = [&](const boost::optional<std::string> & opt) {
      if (opt) readFilepaths_.push_back(datapath + "/" + *opt);
    };
    addFile(params.ocn_file.value());
    addFile(params.ice_file.value());
    addFile(params.fix_file.value());

    // Determine the file to write to (ocean file, or legacy filename)
    if (params.ocn_file.value()) {
      writeFilepath_ = datapath + "/" + *params.ocn_file.value();
    } else if (params.filename.value()) {
      writeFilepath_ = *params.filename.value();
    } else {
      writeFilepath_ = "MOM.res.nc";
    }

    // Backward compatibility: if no named files provided, fall back to filename parameter
    if (readFilepaths_.empty()) {
      readFilepaths_.push_back(writeFilepath_);
    }

    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoMOM6::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                    const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read state");
    oops::Log::trace() << classname() << " read state starting" << std::endl;

    std::vector<std::string> fileVarNames;
    std::vector<double> scalings;
    for (const auto & field : x) {
      const std::string jediName = field.name();
      fileVarNames.push_back(fileionames.has(jediName)
                             ? fileionames.getString(jediName) : "");
      scalings.push_back(fileioscaling.has(jediName)
                         ? fileioscaling.getDouble(jediName) : 0.0);
    }

    readMOM6Netcdf(readFilepaths_, x, fileVarNames, scalings, geom_.getComm());

    oops::Log::trace() << classname() << " read state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoMOM6::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                     const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;

    const auto modelData = geom_.modelData();
    const int ni = modelData.getInt("ni");
    const int nj = modelData.getInt("nj");
    const int nz = modelData.getInt("nz");

    std::vector<std::string> fileVarNames;
    for (const auto & field : x) {
      const std::string jediName = field.name();
      fileVarNames.push_back(fileionames.has(jediName)
                             ? fileionames.getString(jediName) : "");
    }

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

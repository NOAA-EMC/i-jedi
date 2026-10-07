#include "ijedi/Io/cice/IoCICE.h"

#include <fstream>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/cice/UpdateCICENetcdf.h"
#include "ijedi/Io/mom6/ReadMOM6Netcdf.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoCICE> makerIoCICE_("cice");
  // -------------------------------------------------------------------------------------------------
  IoCICE::IoCICE(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()),
        geom_(geom)
  {
    oops::Log::trace() << classname() << " constructor starting" << std::endl;
    const std::string datapath = params.datapath.value();
    filepath_ = datapath + "/" + params.restart_file.value();
    if (params.template_file.value()) {
      templatePath_ = datapath + "/" + *params.template_file.value();
    }
    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoCICE::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
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

    // The MOM6 reader handles the CICE layout: dimensions (ni, nj) and, for fields with several
    // levels, the category dimension
    const auto modelData = geom_.modelData();
    readMOM6Netcdf({filepath_}, x, fileVarNames, scalings,
                   geom_.fields().field("mask3d"), geom_.getFieldMetadata(),
                   modelData.getInt("ni"), modelData.getInt("nj"),
                   modelData.getBool("has_northern_fold", false), geom_.getComm(), {"ncat"});

    oops::Log::trace() << classname() << " read state done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoCICE::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                     const eckit::LocalConfiguration &) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;

    const auto & comm = geom_.getComm();
    if (!templatePath_.empty() && comm.rank() == 0) {
      std::ifstream src(templatePath_, std::ios::binary);
      ASSERT_MSG(src.good(), "IoCICE::write: cannot open template " + templatePath_);
      std::ofstream dst(filepath_, std::ios::binary | std::ios::trunc);
      dst << src.rdbuf();
      ASSERT_MSG(dst.good(), "IoCICE::write: cannot copy the template to " + filepath_);
    }

    std::vector<std::string> fileVarNames;
    for (const auto & field : x) {
      const std::string jediName = field.name();
      fileVarNames.push_back(fileionames.has(jediName)
                             ? fileionames.getString(jediName) : "");
    }
    const auto modelData = geom_.modelData();
    updateCICENetcdf(filepath_, x, fileVarNames, modelData.getInt("ni"), modelData.getInt("nj"),
                     comm);

    oops::Log::trace() << classname() << " write state done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoCICE::print(std::ostream &os) const
  {
    os << classname() << " Io for CICE restarts on the MOM6 grid: " << filepath_;
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

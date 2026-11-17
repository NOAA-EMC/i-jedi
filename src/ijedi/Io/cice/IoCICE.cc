#include "ijedi/Io/cice/IoCICE.h"

#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Geometry/mom6/GeometryMOM6.h"
#include "ijedi/Io/cice/UpdateCICENetcdf.h"
#include "ijedi/Io/mom6/FillMaskedCells.h"
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
    filepath_ = params.datapath.value() + "/" + params.restart_file.value();
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

    // The MOM6 block reader handles the CICE layout: dimensions (ni, nj) and, for fields with
    // several levels, the category dimension.  Each task reads only its MOM6 block.
    const auto * mom6Geom = dynamic_cast<const GeometryMOM6 *>(&geom_.geometryImpl());
    if (mom6Geom == nullptr)
      throw eckit::UserError(classname() + ": the cice io needs a MOM6 geometry", Here());
    readMOM6NetcdfBlocks({filepath_}, x, fileVarNames, scalings, *mom6Geom, {"ncat"});

    // Boundary conditions of masked cells, as for the MOM6 io (the per-category CICE restart
    // fields have none: their masked cells keep the values of the file)
    applyBoundaryConditions(x, geom_.fields().field("mask3d"), geom_.getFieldMetadata());

    oops::Log::trace() << classname() << " read state done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoCICE::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                     const eckit::LocalConfiguration &) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;

    // The restart is updated in place: it must already exist (staged by the workflow)
    const auto & comm = geom_.getComm();
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

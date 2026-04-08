#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Geometry/mpas/GeometryMPAS.h"
#include "ijedi/Io/mpas/IoMPAS.h"
#include "ijedi/Io/mpas/ReadMPASNetcdf.h"
#include "ijedi/Io/mpas/WriteMPASNetcdf.h"

namespace ijedi
{

static IoMaker<IoMPAS> makerIoMPAS_("mpas");

IoMPAS::IoMPAS(const Geometry &geom, const Parameters_ &params)
    : IoBase(geom, params.toConfiguration()),
      geom_(geom),
      datapath_(params.datapath),
      filename_(params.filename.value().value_or(""))
{
  util::Timer timer(classname(), "IoMPAS");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}

IoMPAS::~IoMPAS()
{
  util::Timer timer(classname(), "~IoMPAS");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}

void IoMPAS::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                  const eckit::LocalConfiguration &fileioscaling) const
{
  util::Timer timer(classname(), "read state");
  oops::Log::trace() << classname() << " read state starting" << std::endl;

  const std::string filepath = (!filename_.empty() && filename_[0] == '/')
                                   ? filename_
                                   : datapath_ + "/" + filename_;

  const auto & geomMPAS = dynamic_cast<const GeometryMPAS &>(geom_.geometryImpl());

  std::vector<std::string> fileVarNames;
  std::vector<double> scalings;
  fileVarNames.reserve(x.size());
  scalings.reserve(x.size());

  for (const auto & field : x)
  {
    const std::string jediName = field.name();
    fileVarNames.push_back(fileionames.has(jediName)
                               ? fileionames.getString(jediName)
                               : jediName);
    scalings.push_back(fileioscaling.has(jediName)
                           ? fileioscaling.getDouble(jediName)
                           : 0.0);
  }

  readMPASNetcdf(filepath, x, fileVarNames, scalings, geomMPAS.nCellsGlobal(), geom_.getComm());

  oops::Log::trace() << classname() << " read state done" << std::endl;
}

void IoMPAS::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                   const eckit::LocalConfiguration &fileioscaling) const
{
  util::Timer timer(classname(), "write state");
  oops::Log::trace() << classname() << " write state starting" << std::endl;

  const std::string filepath = (!filename_.empty() && filename_[0] == '/')
                                   ? filename_
                                   : datapath_ + "/" + filename_;

  const auto & geomMPAS = dynamic_cast<const GeometryMPAS &>(geom_.geometryImpl());
  (void)fileioscaling;

  std::vector<std::string> fileVarNames;
  fileVarNames.reserve(x.size());
  for (const auto & field : x)
  {
    const std::string jediName = field.name();
    fileVarNames.push_back(fileionames.has(jediName)
                               ? fileionames.getString(jediName)
                               : jediName);
  }

  writeMPASNetcdf(filepath, x, fileVarNames, geomMPAS.nCellsGlobal(), geomMPAS.nVertLevels(),
                  geom_.getComm());

  oops::Log::trace() << classname() << " write state done" << std::endl;
}

void IoMPAS::print(std::ostream &os) const
{
  os << classname() << " Io for MPAS restart files";
}

}  // namespace ijedi

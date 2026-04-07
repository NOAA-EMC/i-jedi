#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/fv3-history/IoFV3History.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoFV3History> makerIoFV3History_("fv3 history");
  // -------------------------------------------------------------------------------------------------
  IoFV3History::IoFV3History(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration())
  {
    util::Timer timer(classname(), "IoFV3History");
    oops::Log::trace() << classname() << " constructor starting" << std::endl;
    // CALL CONSTRUCTOR
    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  IoFV3History::~IoFV3History()
  {
    util::Timer timer(classname(), "~IoFV3History");
    oops::Log::trace() << classname() << " destructor starting" << std::endl;
    // CALL DESTRUCTOR
    oops::Log::trace() << classname() << " destructor done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoFV3History::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                          const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read state");
    oops::Log::trace() << classname() << " read state starting" << std::endl;
    oops::Log::trace() << classname() << " read state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoFV3History::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                           const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;
    // CALL WRITE
    oops::Log::trace() << classname() << " write state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoFV3History::print(std::ostream &os) const
  {
    os << classname() << " Io for Cube Sphere Histories";
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

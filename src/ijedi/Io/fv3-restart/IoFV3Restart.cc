#include <ostream>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/fv3-restart/IoFV3Restart.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoFV3Restart> makerIoFV3Restart_("fv3 restart");
  // -------------------------------------------------------------------------------------------------
  IoFV3Restart::IoFV3Restart(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration())
  {
    util::Timer timer(classname(), "IoFV3Restart");
    oops::Log::trace() << classname() << " constructor starting" << std::endl;
    // CALL CONSTRUCTOR
    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  IoFV3Restart::~IoFV3Restart()
  {
    util::Timer timer(classname(), "~IoFV3Restart");
    oops::Log::trace() << classname() << " destructor starting" << std::endl;
    // CALL DESTRUCTOR
    oops::Log::trace() << classname() << " destructor done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoFV3Restart::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                          const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read state");
    oops::Log::trace() << classname() << " read state starting" << std::endl;
    oops::Log::trace() << classname() << " read state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoFV3Restart::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                           const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "write state");
    oops::Log::trace() << classname() << " write state starting" << std::endl;
    // CALL WRITE
    oops::Log::trace() << classname() << " write state done" << std::endl;
  }
  // -------------------------------------------------------------------------------------------------
  void IoFV3Restart::print(std::ostream &os) const
  {
    os << classname() << " Io for FMS Restarts";
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

#include <ostream>
#include <string>

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/output/Gmsh.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/gmsh/IoGmsh.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoGmsh> makerIoGmsh_("gmsh");
  // -------------------------------------------------------------------------------------------------
  IoGmsh::IoGmsh(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()),
        geom_(geom),
        filename_(params.filename)
  {
    util::Timer timer(classname(), "IoGmsh");
    oops::Log::trace() << classname() << " constructor starting" << std::endl;
    oops::Log::trace() << classname() << " constructor done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoGmsh::read(atlas::FieldSet &, const eckit::LocalConfiguration &,
                    const eckit::LocalConfiguration &) const
  {
    ABORT("IoGmsh: read is not supported (Gmsh output is write-only)");
  }

  // -------------------------------------------------------------------------------------------------
  void IoGmsh::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &,
                     const eckit::LocalConfiguration &) const
  {
    util::Timer timer(classname(), "write");
    oops::Log::trace() << classname() << " write starting" << std::endl;

    const atlas::functionspace::NodeColumns fspace(geom_.functionSpace());

    eckit::LocalConfiguration gmshConf;
    gmshConf.set("coordinates", "xyz");
    gmshConf.set("ghost", true);
    atlas::output::Gmsh gmsh(filename_, gmshConf);

    gmsh.write(fspace.mesh());
    gmsh.write(x, geom_.functionSpace());

    oops::Log::info() << classname() << " wrote " << filename_ << std::endl;
    oops::Log::trace() << classname() << " write done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoGmsh::print(std::ostream &os) const
  {
    os << classname() << " Gmsh visualization output (" << filename_ << ")";
  }
  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

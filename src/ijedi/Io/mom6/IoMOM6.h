#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/IoBase.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  class IoMOM6Parameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoMOM6Parameters, IoParametersBase)

   public:
    // Path prepended to all files
    oops::Parameter<std::string> datapath{"datapath",
                                          "path to location of files to be read",
                                          "./", this};

    // Ocean background file
    oops::OptionalParameter<std::string> ocn_file{"ocn_file",
                                                  "ocean background file name",
                                                  this};

    // Sea ice file
    oops::OptionalParameter<std::string> ice_file{"ice_file",
                                                  "sea ice file name",
                                                  this};

    // Fix file (decorrelation length scales, distance from coast, etc.)
    oops::OptionalParameter<std::string> fix_file{"fix_file",
                                                  "fix fields file name",
                                                  this};

    // Legacy single filename (kept for backward compatibility)
    oops::OptionalParameter<std::string> filename{"filename",
                                                  "name of the restart or history file (use ocn_file instead)",
                                                  this};
  };

  // -------------------------------------------------------------------------------------------------

  class IoMOM6 : public IoBase, private util::ObjectCounter<IoMOM6>
  {
   public:
    static const std::string classname() { return "ijedi::IoMOM6"; }

    typedef IoMOM6Parameters Parameters_;

    IoMOM6(const Geometry &, const Parameters_ &);
    ~IoMOM6() = default;
    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry & geom_;
    std::vector<std::string> readFilepaths_;  // ordered: ocn, ice, fix
    std::string writeFilepath_;               // ocean file for write
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

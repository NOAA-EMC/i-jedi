#pragma once

#include <ostream>
#include <string>

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
    oops::Parameter<std::string> datapath{"datapath", "path to location of files to be read",
                                          "./", this};

    // Single filename provided
    oops::OptionalParameter<std::string> filename{"filename",
                                                  "name of the restart or history file to be read/written",
                                                  this};
  };

  // -------------------------------------------------------------------------------------------------

  class IoMOM6 : public IoBase, private util::ObjectCounter<IoMOM6>
  {
   public:
    static const std::string classname() { return "ijedi::IoMOM6"; }

    typedef IoMOM6Parameters Parameters_;

    IoMOM6(const Geometry &, const Parameters_ &);
    ~IoMOM6();
    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry & geom_;
    std::string datapath_;
    std::string filename_;
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

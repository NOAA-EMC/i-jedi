#pragma once

#include <ostream>
#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/IoBase.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  class IoCICEParameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoCICEParameters, IoParametersBase)

   public:
    // Path prepended to all files
    oops::Parameter<std::string> datapath{"datapath",
                                          "path to location of files to be read",
                                          "./", this};

    // CICE restart file: read, and updated in place on write
    oops::RequiredParameter<std::string> restart_file{"restart_file",
                                                      "CICE restart file name",
                                                      this};

    // On write, copy this file to the restart file first (the variables not written keep
    // the template's values)
    oops::OptionalParameter<std::string> template_file{"template_file",
                                                       "CICE restart copied to restart_file "
                                                       "before it is updated on write",
                                                       this};
  };

  // -------------------------------------------------------------------------------------------------

  /// Io for CICE restarts on the MOM6 (tripolar) grid of the MOM6 geometry. Fields with several
  /// levels are read from, and written to, variables with an "ncat" category dimension. A write
  /// updates the existing variables of an existing restart (at the geometry's mesh nodes only), so
  /// the variables and grid points that JEDI does not hold keep their values.
  class IoCICE : public IoBase, private util::ObjectCounter<IoCICE>
  {
   public:
    static const std::string classname() { return "ijedi::IoCICE"; }

    typedef IoCICEParameters Parameters_;

    IoCICE(const Geometry &, const Parameters_ &);
    ~IoCICE() = default;
    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry & geom_;
    std::string filepath_;
    std::string templatePath_;
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

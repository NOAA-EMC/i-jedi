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

  class IoGmshParameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoGmshParameters, IoParametersBase)

   public:
    // Output filename (e.g. "state.msh"); Atlas will write per-rank files state.msh.pN
    // and a master state.msh on rank 0 that auto-merges them in the Gmsh GUI.
    oops::RequiredParameter<std::string> filename{"filename",
                                      "name of the Gmsh output file (.msh)",
                                      this};
  };

  // -------------------------------------------------------------------------------------------------

  class IoGmsh : public IoBase, private util::ObjectCounter<IoGmsh>
  {
   public:
    static const std::string classname() { return "ijedi::IoGmsh"; }

    typedef IoGmshParameters Parameters_;

    IoGmsh(const Geometry &, const Parameters_ &);
    ~IoGmsh() = default;
    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry & geom_;
    std::string filename_;
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

#pragma once

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Io/IoBase.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------
  // Parameters used by oops `util::writeFieldSet` and `util::readFieldSet` functions
  // to read/write increments
  class IoUnstructuredParameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoUnstructuredParameters, IoParametersBase)

   public:
    oops::RequiredParameter<std::string> filepath{"filepath",
          "path to the file to be read/written", this};
    oops::Parameter<std::string> extension{"netcdf extension",
          "extension to use for netcdf files", "nc", this};
    oops::Parameter<bool> latStoN{"latitude south to north",
          "whether to order latitudes from south to north", true, this};
    oops::Parameter<bool> oneFilePerTask{"one file per task",
          "whether to create one file per task", false, this};
    oops::Parameter<bool> checkDims{"check dimensions",
          "whether to check dimensions of the file match the geometry", true, this};
  };

  // -------------------------------------------------------------------------------------------------
  class IoUnstructured : public IoBase, private util::ObjectCounter<IoUnstructured>
  {
   public:
    static const std::string classname() { return "ijedi::IoUnstructured"; }

    typedef IoUnstructuredParameters Parameters_;

    IoUnstructured(const Geometry &, const Parameters_ &);
    ~IoUnstructured() = default;
    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry& geometry_;
    const Parameters_ params_;
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

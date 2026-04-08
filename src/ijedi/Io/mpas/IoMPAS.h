#pragma once

#include <ostream>
#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/IoBase.h"

namespace ijedi
{

class IoMPASParameters : public IoParametersBase
{
  OOPS_CONCRETE_PARAMETERS(IoMPASParameters, IoParametersBase)

 public:
  oops::Parameter<std::string> datapath{"datapath", "path to location of files", "./", this};
  oops::OptionalParameter<std::string> filename{"filename",
                                                 "name of the file to read or write",
                                                 this};
};

class IoMPAS : public IoBase, private util::ObjectCounter<IoMPAS>
{
 public:
  static const std::string classname() { return "ijedi::IoMPAS"; }

  typedef IoMPASParameters Parameters_;

  IoMPAS(const Geometry &, const Parameters_ &);
  ~IoMPAS();

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

}  // namespace ijedi

#include "ijedi/Io/unstructured/IoUnstructured.h"

#include <vector>

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

static IoMaker<IoUnstructured> makerIoUnstructured_("unstructured");

// -------------------------------------------------------------------------------------------------

IoUnstructured::IoUnstructured(const Geometry &geom, const Parameters_ &params)
    : IoBase(geom, params.toConfiguration()), geometry_(geom), params_(params)
{}

// -------------------------------------------------------------------------------------------------

void IoUnstructured::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                          const eckit::LocalConfiguration &fileioscaling) const {
  util::Timer timer(classname(), "read state");
  oops::Log::trace() << classname() << " read state starting" << std::endl;

  // TODO(AS): clean this up once it's more clear where variables vs levels vs metadata
  // are evolving
  oops::Variables varsToRead;
  for (const auto & field : x) {
    const std::string & currentVar = field.name();
    if (fileionames.has(currentVar)) {
      varsToRead.push_back(fileionames.getString(currentVar));
    } else {
      varsToRead.push_back(currentVar);
    }
  }
  std::vector<size_t> varSizes = geometry_.variableSizes(varsToRead);
  x.clear();
  oops::Log::info() << "Reading variables: " << varsToRead << std::endl;
  util::readFieldSet(geometry_.comm(), geometry_.functionSpace(), varSizes,
                     varsToRead.variables(), params_.toConfiguration(), x);
  oops::Log::trace() << classname() << " read state done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void IoUnstructured::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                           const eckit::LocalConfiguration &fileioscaling) const {
  util::Timer timer(classname(), "write state");
  oops::Log::trace() << classname() << " write state starting" << std::endl;
  util::writeFieldSet(geometry_.comm(), params_.toConfiguration(), x);
  oops::Log::trace() << classname() << " write state done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void IoUnstructured::print(std::ostream &os) const {
  os << classname() << " oops utils implementation of read/write for unstructured grids";
}

}  // namespace ijedi

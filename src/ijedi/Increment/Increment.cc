#include "ijedi/Increment/Increment.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>

#include "atlas/field.h"
#include "eckit/config/Configuration.h"
#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Utilities/PrintHelper.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/for_each.h"

namespace ijedi {

  // -----------------------------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom, const oops::Variables & vars,
                       const util::DateTime & time)
      : mist::base::Increment(geom, vars, time), geom_(geom) {}

  // -----------------------------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom, const Increment & other, const bool ad)
      : mist::base::Increment(geom, other, ad), geom_(geom) {}

  // -----------------------------------------------------------------------------------------------

  Increment::Increment(const Increment & other, const bool copy)
      : mist::base::Increment(other, copy), geom_(other.geom_) {}

  // -----------------------------------------------------------------------------------------------

  Increment::~Increment() = default;

  // -----------------------------------------------------------------------------------------------

  Increment & Increment::operator=(const Increment & rhs) {
    mist::base::Increment::operator=(rhs);
    return *this;
  }

  // -----------------------------------------------------------------------------------------------

  void Increment::read(const eckit::Configuration & config) {
    oops::Log::trace() << "ijedi::Increment::read starting" << std::endl;

    // Create a Parameters object
    IncrementParameters params;
    params.deserialize(config);

    // Check that there are IO parameters
    if (params.io.value() == boost::none ||
        params.io.value()->ioParameters.value() == nullptr)
    {
      throw eckit::BadParameter("ijedi::Increment::read: No IO parameters provided", Here());
    }

    // Get the polymorphic IO parameters
    const IoParametersBase &ioParams = *params.io.value()->ioParameters.value();

    // Create the IO object to use
    // ---------------------------
    std::unique_ptr<IoBase> io(IoFactory::create(geom_, ioParams));

    // Call read method of child
    // -------------------------
    io->readBase(this->fieldSet());

    oops::Log::trace() << "ijedi::Increment::read done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  void Increment::write(const eckit::Configuration & config) const {
    oops::Log::trace() << "ijedi::Increment::write starting" << std::endl;

    // Create a Parameters object
    IncrementParameters params;
    params.deserialize(config);

    // Check that there are IO parameters
    if (params.io.value() == boost::none ||
        params.io.value()->ioParameters.value() == nullptr)
    {
      throw eckit::BadParameter("ijedi::Increment::read: No IO parameters provided", Here());
    }

    // Get the polymorphic IO parameters
    const IoParametersBase &ioParams = *params.io.value()->ioParameters.value();

    // Create the IO object to use
    // ---------------------------
    std::unique_ptr<IoBase> io(IoFactory::create(geom_, ioParams));

    // Call write method of child
    // --------------------------
    io->writeBase(this->fieldSet());

    oops::Log::trace() << "ijedi::Increment::write done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  void Increment::dirac(const eckit::Configuration & config) {
    oops::Log::trace() << "ijedi::Increment::dirac starting" << std::endl;

    // Create a Parameters object
    DiracParameters params;
    params.deserialize(config);

    // Extract parameters
    const std::vector<std::string> diracFlds = params.diracFlds;
    const std::vector<int> diracProc = params.diracProc;
    const std::vector<int> diracHorx = params.diracHorx;
    const std::vector<int> diracVert = params.diracVert;

    // Assert that all vectors are the same lenght
    size_t nDiracs = diracFlds.size();
    ASSERT_MSG(diracProc.size() == nDiracs, "Dirac parameters diracProc incorrect length");
    ASSERT_MSG(diracHorx.size() == nDiracs, "Dirac parameters diracHorx incorrect length");
    ASSERT_MSG(diracVert.size() == nDiracs, "Dirac parameters diracVert incorrect length");

    // Set diracs
    for (size_t i = 0; i < diracFlds.size(); ++i) {
      const std::string & varName = diracFlds[i];
      const int proc = diracProc[i];
      const int horx = diracHorx[i];
      const int vert = diracVert[i];

      if (this->geom_.comm().rank() == proc) {
        auto field = this->fieldSet().field(varName);
        auto view = atlas::array::make_view<double, 2>(field);
        view(horx, vert) = 1.0;
      }
    }
    oops::Log::trace() << "ijedi::Increment::dirac done" << std::endl;
  }

  // -----------------------------------------------------------------------------------------------

  void Increment::print(std::ostream & os) const {
    os << std::endl
       << "  Valid time: " << this->validTime() << std::endl
       << ", nFields = " << this->variables().size();

    const auto & comm = geom_.comm();
    const auto & fs   = this->fieldSet();
    for (const auto & var : this->variables()) {
      const atlas::Field & field            = fs.field(var.name());
      const auto[globalMin, globalMax, rms] = fieldMinMaxRMS(comm, field);
      os << std::endl
         << var.name() << " : " << std::scientific << std::setprecision(16) << "Min=" << globalMin
         << ", Max=" << globalMax << ", RMS=" << rms;
    }
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi

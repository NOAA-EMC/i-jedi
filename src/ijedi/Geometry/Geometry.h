#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "eckit/mpi/Comm.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "mist/base/Geometry.h"

#include "ijedi/Geometry/base/GeometryBase.h"

#include "ijedi/FieldMetadata/FieldsMetadata.h"

// Forward declarations
namespace eckit
{
  class Configuration;
}

namespace oops
{
  class Variables;
}

namespace ijedi
{

  // -----------------------------------------------------------------------------
  // Geometry handles geometry.

  class Geometry : public mist::Geometry, private util::ObjectCounter<Geometry>
  {
   public:
    static const std::string classname() { return "ijedi::Geometry"; }

    Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
    ~Geometry();

    // This might need to change with mist.
    const int &numLevels() const { return numberLevels_; }

    // Function to access field metadata
    const FieldsMetadata &getFieldMetadata() const { return *fieldsMeta_; }

    // Add model-specific Vader ingredient fields (masks, area fractions, ...)
    // to fset, sourced from the built geometry fields. Coordinates are handled
    // generically by the caller (see addVaderGeometryIngredients).
    void addModelVaderIngredients(atlas::FieldSet &fset) const
      { geometryImpl_->addVaderIngredients(fields(), fset, numberLevels_); }

   private:
    Geometry &operator=(const Geometry &);
    void print(std::ostream &) const;
    std::shared_ptr<FieldsMetadata> fieldsMeta_;
    std::shared_ptr<GeometryBase> geometryImpl_;
    int numberLevels_;
    std::string type_;
  };
  // -----------------------------------------------------------------------------

}  // namespace ijedi

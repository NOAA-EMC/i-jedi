#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class GeometryBase
  {
   public:
    virtual ~GeometryBase() = default;

    static std::shared_ptr<GeometryBase> create(const eckit::Configuration &,
                                                const eckit::mpi::Comm &,
                                                eckit::LocalConfiguration &,
                                                atlas::FunctionSpace &,
                                                atlas::FieldSet &,
                                                bool &, int &);
    virtual void print(std::ostream &) const = 0;
    virtual std::vector<double> verticalCoord(std::string &) const = 0;

    /// \brief Add model-specific Vader ingredient fields (e.g. masks, area
    ///        fractions) to \p fset, sourced from the already-built geometry
    ///        fields \p geomFields. Coordinates (latitude/longitude) are handled
    ///        generically by the caller and are not the responsibility of models.
    ///        Default: no-op (models without ocean-style ingredients need nothing).
    virtual void addVaderIngredients(const atlas::FieldSet & geomFields,
                                     atlas::FieldSet & fset, int nlevels) const {}
  };

}  // namespace ijedi

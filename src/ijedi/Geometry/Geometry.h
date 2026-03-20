// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

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

  // ---------------------------------------------------------------------------
  // Geometry handles geometry.

  class Geometry : public util::Printable,
                   private util::ObjectCounter<Geometry>
  {
   public:
    static const std::string classname() { return "ijedi::Geometry"; }

    explicit Geometry(const eckit::Configuration &,
                      const eckit::mpi::Comm &);
    Geometry(const Geometry &);
    ~Geometry();

    bool levelsAreTopDown() const { return true; }
    std::vector<double> verticalCoord(std::string &) const;
    std::vector<size_t> variableSizes(const oops::Variables &) const;

    const eckit::mpi::Comm &getComm() const { return comm_; }

    const atlas::FunctionSpace &functionSpace() const
      { return geometryImpl_->functionSpace(); }
    const atlas::FieldSet &fields() const
      { return geometryImpl_->fields(); }
    atlas::FunctionSpace &functionSpace()
      { return geometryImpl_->functionSpace(); }
    atlas::FieldSet &fields()
      { return geometryImpl_->fields(); }
    const int &numLevels() const { return geometryImpl_->numLevels(); }

    // Access to grid-specific parameters
    eckit::LocalConfiguration gridSpecific() const
      { return geometryImpl_->gridSpecific(); }

   private:
    Geometry &operator=(const Geometry &);
    void print(std::ostream &) const;
    const eckit::mpi::Comm &comm_;
    std::shared_ptr<FieldsMetadata> fieldsMeta_;
    std::shared_ptr<GeometryBase> geometryImpl_;
  };
  // ---------------------------------------------------------------------------

}  // namespace ijedi

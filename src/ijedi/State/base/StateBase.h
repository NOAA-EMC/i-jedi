// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <memory>
#include <ostream>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class Geometry;

  // ---------------------------------------------------------------------------
  // Abstract base class for model-specific State implementations.
  //
  // Mirrors GeometryBase: subclasses populate functionSpace_ and numLevels_
  // in their constructors by calling one of the protected init* helpers,
  // which set the common data members and allocate fields.  All generic
  // Atlas-based operations (zero, accumul, norm, serialize, etc.) are
  // implemented once here.  Only print, read, and write are pure virtual.
  //
  // A model-specific subclass constructor is typically two lines:
  //   : StateBase(geom.getComm()), geom_(geom)
  //   { initFromGeomVarsTime(geom, vars, time); }
  // ---------------------------------------------------------------------------
  class StateBase
  {
   public:
    virtual ~StateBase() = default;

    // Factory methods: dispatch on geom.gridSpecific().getString("grid_type")
    static std::shared_ptr<StateBase> create(const Geometry &,
                                             const oops::Variables &,
                                             const util::DateTime &);
    static std::shared_ptr<StateBase> create(const Geometry &,
                                             const eckit::Configuration &);
    // change geometry:
    static std::shared_ptr<StateBase> create(const Geometry &,
                                             const StateBase &);
    // change variables:
    static std::shared_ptr<StateBase> create(const oops::Variables &,
                                             const StateBase &);

    // Non-virtual accessors
    const oops::Variables  & variables() const { return vars_; }
    const util::DateTime   & validTime() const { return time_; }
    void updateTime(const util::Duration & dt) { time_ += dt; }

    const atlas::FieldSet  & fields() const { return fields_; }
    atlas::FieldSet        & fields()       { return fields_; }

    // Copy data (vars, time, field values) from another state of the same
    // geometry.  Used by State::operator= — does not touch functionSpace_,
    // numLevels_, or comm_.
    void copyDataFrom(const StateBase &);

    // Generic Atlas-based operations — implemented once in StateBase
    void zero();
    void accumul(const double &, const StateBase &);
    double norm() const;
    size_t serialSize() const;
    void serialize(std::vector<double> &) const;
    void deserialize(const std::vector<double> &, size_t &);
    void toFieldSet(atlas::FieldSet &) const;
    void fromFieldSet(const atlas::FieldSet &);

    // Default implementation prints per-variable min/max/rms stats.
    // Subclasses may override to add model-specific information.
    virtual void print(std::ostream &) const;

    // Pure virtual — file format is model-specific
    virtual void read(const eckit::Configuration &) = 0;
    virtual void write(const eckit::Configuration &) const = 0;

   protected:
    // The only constructor — subclass must pass the communicator
    explicit StateBase(const eckit::mpi::Comm & comm) : comm_(comm) {}

    // Init helpers — called from subclass constructors to set common data
    // and allocate fields.  Each sets functionSpace_, numLevels_, vars_,
    // time_, then calls allocateFields() (with optional field copy for
    // change-vars).
    void initFromGeomVarsTime(const Geometry &, const oops::Variables &,
                              const util::DateTime &);
    void initFromConf(const Geometry &, const eckit::Configuration &);
    void initFromGeomState(const Geometry &, const StateBase &);
    void initFromVarsState(const Geometry &, const oops::Variables &,
                           const StateBase &);

    // Data members set via init helpers — mirrors GeometryBase layout
    atlas::FunctionSpace   functionSpace_;
    atlas::FieldSet        fields_;
    oops::Variables        vars_;
    util::DateTime         time_;
    int                    numLevels_ = 0;
    const eckit::mpi::Comm & comm_;

   private:
    void allocateFields();
  };

}  // namespace ijedi

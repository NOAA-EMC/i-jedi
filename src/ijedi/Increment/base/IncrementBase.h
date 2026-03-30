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
  class StateBase;

  // ---------------------------------------------------------------------------
  // Abstract base class for model-specific Increment implementations.
  //
  // Mirrors StateBase: subclasses populate functionSpace_ and numLevels_
  // in their constructors by calling one of the protected init* helpers.
  // All generic Atlas-based linear-algebra operations are implemented once
  // here.  Only print, read, and write are pure virtual.
  //
  // A model-specific subclass constructor is typically two lines:
  //   : IncrementBase(geom.getComm()), geom_(geom)
  //   { initFromGeomVarsTime(geom, vars, time); }
  // ---------------------------------------------------------------------------
  class IncrementBase
  {
   public:
    virtual ~IncrementBase() = default;

    // Factory methods: dispatch on geom.gridSpecific().getString("grid_type")
    static std::shared_ptr<IncrementBase> create(const Geometry &,
                                                 const oops::Variables &,
                                                 const util::DateTime &);
    // change resolution:
    static std::shared_ptr<IncrementBase> create(const Geometry &,
                                                 const IncrementBase &);
    // copy (with optional zero):
    static std::shared_ptr<IncrementBase> create(const IncrementBase &,
                                                 bool copy);

    // Non-virtual accessors
    const oops::Variables  & variables() const { return vars_; }
    const util::DateTime   & validTime() const { return time_; }
    void updateTime(const util::Duration & dt) { time_ += dt; }

    const atlas::FieldSet  & fields() const { return fields_; }
    atlas::FieldSet        & fields()       { return fields_; }

    // Copy data from another increment of the same geometry (used by operator=)
    void copyDataFrom(const IncrementBase &);

    // -----------------------------------------------------------------------
    // Generic Atlas-based operations — implemented once in IncrementBase
    // -----------------------------------------------------------------------

    // Inherited from State pattern
    void zero();
    void zero(const util::DateTime & dt) { zero(); time_ = dt; }
    double norm() const;
    size_t serialSize() const;
    void serialize(std::vector<double> &) const;
    void deserialize(const std::vector<double> &, size_t &);
    void toFieldSet(atlas::FieldSet &) const;
    void fromFieldSet(const atlas::FieldSet &);

    // Increment-specific
    void ones();
    void sqrt();
    void random();
    void dirac(const eckit::Configuration &);

    void diff(const StateBase &, const StateBase &);
    void accumul(const double &, const StateBase &);

    IncrementBase & operator+=(const IncrementBase &);
    IncrementBase & operator-=(const IncrementBase &);
    IncrementBase & operator*=(const double &);
    void axpy(const double &, const IncrementBase &);
    double dot_product_with(const IncrementBase &) const;
    void schur_product_with(const IncrementBase &);

    // Default print: per-variable min/max/rms
    virtual void print(std::ostream &) const;

    // Pure virtual — file format is model-specific
    virtual void read(const eckit::Configuration &) = 0;
    virtual void write(const eckit::Configuration &) const = 0;

   protected:
    explicit IncrementBase(const eckit::mpi::Comm & comm) : comm_(comm) {}

    // Init helpers — called from subclass constructors
    void initFromGeomVarsTime(const Geometry &, const oops::Variables &,
                              const util::DateTime &);
    void initFromGeomIncrement(const Geometry &, const IncrementBase &);
    void initCopy(const IncrementBase &, bool copy);

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

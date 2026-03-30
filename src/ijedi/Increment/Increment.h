// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "ijedi/Increment/base/IncrementBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class Geometry;
  class State;

  // ---------------------------------------------------------------------------
  // Increment — OOPS-facing wrapper.
  // Mirrors State: holds a shared_ptr<IncrementBase> and delegates everything.
  // ---------------------------------------------------------------------------
  class Increment : public util::Printable,
                    public util::Serializable,
                    private util::ObjectCounter<Increment>
  {
   public:
    static const std::string classname() { return "ijedi::Increment"; }

    // OOPS-required constructors
    Increment(const Geometry &, const oops::Variables &, const util::DateTime &);
    Increment(const Geometry &, const Increment &);     // change geometry
    Increment(const Increment &, bool copy = true);     // copy / nocopy
    ~Increment();

    Increment & operator=(const Increment &);

    // Time
    const util::DateTime & validTime() const;
    void updateTime(const util::Duration &);

    // Increment operations
    void diff(const State &, const State &);
    void zero();
    void zero(const util::DateTime &);
    void ones();
    void sqrt();
    void dirac(const eckit::Configuration &);
    void random();
    void accumul(const double &, const State &);

    // Linear algebra
    Increment & operator+=(const Increment &);
    Increment & operator-=(const Increment &);
    Increment & operator*=(const double &);
    void axpy(const double &, const Increment &, const bool check = true);
    double dot_product_with(const Increment &) const;
    void schur_product_with(const Increment &);

    // I/O
    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &) const;
    double norm() const;

    // Serialization
    size_t serialSize() const override;
    void serialize(std::vector<double> &) const override;
    void deserialize(const std::vector<double> &, size_t &) override;

    // Atlas FieldSet interface
    void toFieldSet(atlas::FieldSet &) const;
    void fromFieldSet(const atlas::FieldSet &);

    // Variables
    const oops::Variables & variables() const;

    // Internal accessor for IncrementBase (used by other Increment methods)
    const IncrementBase & impl() const { return *impl_; }
    IncrementBase & impl()             { return *impl_; }

   private:
    void print(std::ostream &) const override;
    std::shared_ptr<IncrementBase> impl_;
  };

}  // namespace ijedi

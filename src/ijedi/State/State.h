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

#include "ijedi/State/base/StateBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class Geometry;
  class Increment;

  // ---------------------------------------------------------------------------
  // State — OOPS-facing wrapper.
  // Mirrors Geometry: holds a shared_ptr<StateBase> and delegates everything.
  // ---------------------------------------------------------------------------
  class State : public util::Printable,
                private util::ObjectCounter<State>
  {
   public:
    static const std::string classname() { return "ijedi::State"; }

    // OOPS-required constructors
    State(const Geometry &, const oops::Variables &, const util::DateTime &);
    State(const Geometry &, const eckit::Configuration &);
    State(const Geometry &, const State &);          // change geometry
    State(const oops::Variables &, const State &);   // change variables
    State(const State &);
    ~State();
    State & operator=(const State &);

    // Time
    const util::DateTime & validTime() const;
    void updateTime(const util::Duration &);

    // State operations
    void zero();
    void accumul(const double &, const State &);
    double norm() const;
    State & operator+=(const Increment &);

    // I/O
    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &) const;

    // Serialization (used in 4DEnVar, weak-constraint 4DVar)
    size_t serialSize() const;
    void serialize(std::vector<double> &) const;
    void deserialize(const std::vector<double> &, size_t &);

    // Atlas FieldSet interface
    void toFieldSet(atlas::FieldSet &) const;
    void fromFieldSet(const atlas::FieldSet &);

    // Variables
    const oops::Variables & variables() const;

    // Internal accessor for StateBase (used by Increment::diff and accumul)
    const StateBase & impl() const { return *stateImpl_; }

   private:
    void print(std::ostream &) const override;
    std::shared_ptr<StateBase> stateImpl_;
  };

}  // namespace ijedi

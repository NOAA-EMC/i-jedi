/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/IoBase.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  class IoSQGParameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoSQGParameters, IoParametersBase)
   public:
    oops::RequiredParameter<std::string> filename{"filename",
      "SQG netCDF file (sqgturb-cpp layout: pv(t,z,y,x) in m/s). When writing, \"%{valid}\" "
      "in the name is replaced by the valid time (e.g. 20150101T060000Z)", this};
    oops::Parameter<int> timeIndex{"time index",
      "record of the file to read (negative counts from the end; -1 is the last)", -1, this};
    oops::Parameter<float> scaleToK{"scale to K",
      "f*theta0/g used when writing (reading uses the file's scale_to_K attribute)",
      1.e-4f * 300.f / 9.8f, this};
  };

  // -------------------------------------------------------------------------------------------------
  /// Reads/writes the SQG model's netCDF file: boundary PV pv(t,z,y,x) in m/s (z=0 lower
  /// boundary, z=1 lid), converted to/from the state variable "air_potential_temperature" (K)
  /// with the file's "scale_to_K" attribute (f*theta0/g; for older files, computed from the
  /// global attributes f, theta0 and g). The geometry must be the SQG equatorial patch
  /// (see Model/sqg/SQGGrid.h). Writing needs a single MPI task.
  class IoSQG : public IoBase, private util::ObjectCounter<IoSQG>
  {
   public:
    static const std::string classname() { return "ijedi::IoSQG"; }
    typedef IoSQGParameters Parameters_;

    IoSQG(const Geometry &, const Parameters_ &);
    ~IoSQG() = default;

    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry & geom_;
    std::string filename_;
    int timeIndex_;
    float scaleToK_;
  };

  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

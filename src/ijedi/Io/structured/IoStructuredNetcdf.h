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
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Io/IoBase.h"

namespace ijedi
{

  class Geometry;

  // -------------------------------------------------------------------------------------------------

  class IoStructuredNetcdfParameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoStructuredNetcdfParameters, IoParametersBase)

   public:
    oops::RequiredParameter<std::string> filepath{"filepath",
          "path to the netCDF file to be read or written", this};
    oops::Parameter<int> timeIndex{"time index",
          "record to read when the file has a time dimension", 0, this};
    oops::OptionalParameter<std::string> datetime{"datetime",
          "valid time to stamp on the file when writing; without it no time dimension is "
          "written, which is preferable to inventing a wrong timestamp", this};
    oops::Parameter<bool> latStoN{"latitude south to north",
          "latitude ordering to use when writing. The default matches ARCO-ERA5 and "
          "WeatherBench-2. On read the ordering is detected from the file, not from this.",
          true, this};
    oops::Parameter<std::string> levelUnits{"level units",
          "units for the level coordinate when writing: 'hPa' (as in ERA5) or 'Pa'. On read "
          "the units are taken from the file.", "hPa", this};
    oops::Parameter<std::string> precision{"precision",
          "on-disk type for the data variables: 'float' (as in ERA5) or 'double'",
          "double", this};
  };

  // -------------------------------------------------------------------------------------------------

  /// \brief Read and write structured latitude-longitude netCDF files laid out the way the
  ///        AI weather models expect their inputs.
  ///
  /// The layout is the ARCO-ERA5 / WeatherBench-2 convention that NeuralGCM, GraphCast, Pangu
  /// and Aurora all consume: one netCDF variable per field, dimensioned
  /// (time, level, latitude, longitude) for upper-air fields and (time, latitude, longitude)
  /// for surface fields, with CF coordinate variables and SI units. Nothing here is specific
  /// to one model - which variables are exchanged, and under which names, is decided by the
  /// caller through IoBase's "field io names" and "field io scaling" maps.
  ///
  /// Parallelism. Reads are independent and fully distributed: each task opens the file
  /// read-only and reads just the bounding box of the grid it owns, so no task ever holds the
  /// global field and the file access itself involves no collectives. The one collective on
  /// the read path is the halo exchange at the end, which fills the halo points that no task
  /// reads for itself. Writes gather to the root
  /// task, which is the simpler and more portable choice and costs little at the resolutions
  /// these models run at (a 128x64x37 field is 2.4 MB); revisit with parallel netCDF only if
  /// output ever becomes the bottleneck.
  ///
  /// The reader validates the file's grid against the geometry - longitudes, latitudes
  /// (allowing for the two conventional orderings) and pressure levels (allowing hPa or Pa) -
  /// and refuses to read a file that does not match. A silently transposed or flipped field is
  /// the failure mode that costs the most to find later, so it is turned into a startup error.
  class IoStructuredNetcdf : public IoBase, private util::ObjectCounter<IoStructuredNetcdf>
  {
   public:
    static const std::string classname() { return "ijedi::IoStructuredNetcdf"; }

    typedef IoStructuredNetcdfParameters Parameters_;

    IoStructuredNetcdf(const Geometry &, const Parameters_ &);
    ~IoStructuredNetcdf() = default;

    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry & geometry_;
    const Parameters_ params_;
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

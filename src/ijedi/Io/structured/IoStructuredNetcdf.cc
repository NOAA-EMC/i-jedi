/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Io/structured/IoStructuredNetcdf.h"

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/option.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/structured/StructuredGridLayout.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  static IoMaker<IoStructuredNetcdf> makerIoStructuredNetcdf_("structured netcdf");

  // -------------------------------------------------------------------------------------------------

  namespace
  {

    // Coordinate names as they appear in ARCO-ERA5 and WeatherBench-2, plus the short spellings
    // that ECMWF tooling emits.
    const std::vector<std::string> kLonNames = {"longitude", "lon", "x"};
    const std::vector<std::string> kLatNames = {"latitude", "lat", "y"};
    const std::vector<std::string> kLevelNames = {"level", "pressure_level", "plev", "isobaric"};
    const std::vector<std::string> kTimeNames = {"time", "valid_time"};

    // Grids are compared in degrees, so this is generous next to the spacing of even a 0.25
    // degree grid while still catching a half-cell offset or a flipped axis.
    constexpr double kCoordTolerance = 1.0e-4;
    // Pressure levels are compared in Pa; 0.5 Pa separates any two ERA5 levels comfortably.
    constexpr double kLevelTolerance = 0.5;

    void ncCheck(int status, const std::string & what)
    {
      if (status != NC_NOERR)
      {
        throw eckit::Exception("IoStructuredNetcdf: " + what + ": " + nc_strerror(status), Here());
      }
    }

    /// Find the first of \p candidates that exists as a dimension, returning its length.
    /// Returns false when none is present, which is not always an error (time is optional).
    bool findDim(int ncid, const std::vector<std::string> & candidates,
                 std::string * name, int * dimid, size_t * len)
    {
      for (const std::string & candidate : candidates)
      {
        if (nc_inq_dimid(ncid, candidate.c_str(), dimid) == NC_NOERR)
        {
          ncCheck(nc_inq_dimlen(ncid, *dimid, len), "reading length of dimension " + candidate);
          *name = candidate;
          return true;
        }
      }
      return false;
    }

    std::string joinNames(const std::vector<std::string> & names)
    {
      std::string joined;
      for (const std::string & name : names)
      {
        if (!joined.empty()) joined += ", ";
        joined += name;
      }
      return joined;
    }

    std::vector<double> readCoordVar(int ncid, const std::string & name, size_t len)
    {
      int varid;
      ncCheck(nc_inq_varid(ncid, name.c_str(), &varid),
              "the file has a '" + name + "' dimension but no matching coordinate variable, so "
              "its grid cannot be checked against the geometry");
      std::vector<double> values(len);
      ncCheck(nc_get_var_double(ncid, varid, values.data()), "reading coordinate " + name);
      return values;
    }

    /// Read a variable's units attribute, empty when it has none.
    std::string readUnits(int ncid, const std::string & varName)
    {
      int varid;
      if (nc_inq_varid(ncid, varName.c_str(), &varid) != NC_NOERR) return "";
      size_t len = 0;
      if (nc_inq_attlen(ncid, varid, "units", &len) != NC_NOERR) return "";
      std::string units(len, '\0');
      if (nc_get_att_text(ncid, varid, "units", &units[0]) != NC_NOERR) return "";
      units.resize(strnlen(units.c_str(), len));
      return units;
    }

    /// Convert a level coordinate to Pa. ERA5 files label levels in hPa; some tools write Pa.
    /// With no units attribute, fall back on magnitude: no atmospheric level sits above
    /// 2000 hPa, and none below 2000 Pa is expressible in whole hPa in these files.
    std::vector<double> levelsToPa(const std::vector<double> & raw, const std::string & units,
                                   std::string * resolvedUnits)
    {
      const bool isHectoPascal =
          units.find("hPa") != std::string::npos ||
          units.find("millibar") != std::string::npos ||
          units.find("mbar") != std::string::npos ||
          (units.empty() && !raw.empty() &&
           *std::max_element(raw.begin(), raw.end()) < 2000.0);

      *resolvedUnits = isHectoPascal ? "hPa" : "Pa";
      std::vector<double> pa(raw.size());
      for (size_t k = 0; k < raw.size(); ++k)
      {
        pa[k] = isHectoPascal ? raw[k] * 100.0 : raw[k];
      }
      return pa;
    }

    /// How a file's axes line up with the geometry's.
    struct FileMapping
    {
      bool hasTime = false;
      size_t nTime = 0;
      size_t nLevelFile = 0;
      bool latFlipped = false;          ///< file latitudes run opposite to the atlas grid
      std::vector<int> levelIndex;      ///< geometry level -> index along the file level axis
      std::string lonDim, latDim, levelDim, timeDim;
    };

    /// Check that the file describes the same grid as the geometry, and work out the index
    /// mapping. Everything that can go silently wrong - a transposed field, a flipped
    /// hemisphere, levels in the wrong order or the wrong units - is caught here rather than
    /// showing up as a strange forecast days later.
    FileMapping matchFile(int ncid, const StructuredGridLayout & layout,
                          const std::vector<double> & geomLevelsPa, int requestedTimeIndex,
                          const std::string & filepath)
    {
      FileMapping mapping;
      int dimid = 0;
      size_t len = 0;

      // --- longitude ---
      if (!findDim(ncid, kLonNames, &mapping.lonDim, &dimid, &len))
      {
        throw eckit::BadValue("IoStructuredNetcdf: " + filepath + " has no longitude dimension "
                              "(tried: " + joinNames(kLonNames) + ")", Here());
      }
      if (static_cast<int>(len) != layout.nx)
      {
        std::stringstream msg;
        msg << "IoStructuredNetcdf: " << filepath << " has " << len << " longitudes but the "
            << "geometry has " << layout.nx;
        throw eckit::BadValue(msg.str(), Here());
      }
      const std::vector<double> fileLons = readCoordVar(ncid, mapping.lonDim, len);
      for (int i = 0; i < layout.nx; ++i)
      {
        if (std::abs(normaliseLongitude(fileLons[i]) - layout.lons[i]) > kCoordTolerance)
        {
          std::stringstream msg;
          msg << "IoStructuredNetcdf: " << filepath << " longitude " << i << " is "
              << fileLons[i] << " but the geometry expects " << layout.lons[i]
              << ". The file must be on the geometry's grid; regrid it first.";
          throw eckit::BadValue(msg.str(), Here());
        }
      }

      // --- latitude, in either of the two conventional orderings ---
      if (!findDim(ncid, kLatNames, &mapping.latDim, &dimid, &len))
      {
        throw eckit::BadValue("IoStructuredNetcdf: " + filepath + " has no latitude dimension "
                              "(tried: " + joinNames(kLatNames) + ")", Here());
      }
      if (static_cast<int>(len) != layout.ny)
      {
        std::stringstream msg;
        msg << "IoStructuredNetcdf: " << filepath << " has " << len << " latitudes but the "
            << "geometry has " << layout.ny;
        throw eckit::BadValue(msg.str(), Here());
      }
      const std::vector<double> fileLats = readCoordVar(ncid, mapping.latDim, len);

      bool sameOrder = true;
      bool reversedOrder = true;
      for (int j = 0; j < layout.ny; ++j)
      {
        if (std::abs(fileLats[j] - layout.lats[j]) > kCoordTolerance) sameOrder = false;
        if (std::abs(fileLats[layout.ny - 1 - j] - layout.lats[j]) > kCoordTolerance)
        {
          reversedOrder = false;
        }
      }
      if (!sameOrder && !reversedOrder)
      {
        std::stringstream msg;
        msg << "IoStructuredNetcdf: the latitudes in " << filepath << " match the geometry in "
            << "neither order (file runs " << fileLats.front() << " to " << fileLats.back()
            << ", geometry runs " << layout.lats.front() << " to " << layout.lats.back()
            << "). The file must be on the geometry's grid; regrid it first.";
        throw eckit::BadValue(msg.str(), Here());
      }
      mapping.latFlipped = !sameOrder;

      // --- levels ---
      if (findDim(ncid, kLevelNames, &mapping.levelDim, &dimid, &mapping.nLevelFile))
      {
        std::string resolvedUnits;
        const std::vector<double> fileLevelsPa =
            levelsToPa(readCoordVar(ncid, mapping.levelDim, mapping.nLevelFile),
                       readUnits(ncid, mapping.levelDim), &resolvedUnits);

        // Match by value rather than by position, so a file ordered surface-to-top reads
        // correctly against a geometry ordered top-to-surface.
        mapping.levelIndex.assign(geomLevelsPa.size(), -1);
        for (size_t k = 0; k < geomLevelsPa.size(); ++k)
        {
          for (size_t kf = 0; kf < fileLevelsPa.size(); ++kf)
          {
            if (std::abs(fileLevelsPa[kf] - geomLevelsPa[k]) <= kLevelTolerance)
            {
              mapping.levelIndex[k] = static_cast<int>(kf);
              break;
            }
          }
          if (mapping.levelIndex[k] < 0)
          {
            std::stringstream msg;
            msg << "IoStructuredNetcdf: " << filepath << " (levels in " << resolvedUnits
                << ") has no level at " << geomLevelsPa[k] << " Pa, which the geometry needs "
                << "for level " << k;
            throw eckit::BadValue(msg.str(), Here());
          }
        }
      }

      // --- time, optional ---
      if (findDim(ncid, kTimeNames, &mapping.timeDim, &dimid, &mapping.nTime))
      {
        mapping.hasTime = true;
        if (requestedTimeIndex < 0 || static_cast<size_t>(requestedTimeIndex) >= mapping.nTime)
        {
          std::stringstream msg;
          msg << "IoStructuredNetcdf: 'time index' " << requestedTimeIndex << " is out of range "
              << "for " << filepath << ", which has " << mapping.nTime << " time records";
          throw eckit::BadValue(msg.str(), Here());
        }
      }

      return mapping;
    }

    /// The geometry's reference pressure column, in Pa. Empty when the geometry carries no
    /// vertical coordinate, which is legal for a surface-only file.
    std::vector<double> geometryLevelsPa(const Geometry & geom)
    {
      const eckit::LocalConfiguration modelData = geom.modelData();
      if (!modelData.has("vertical_coordinate_reference_pressure")) return {};
      return modelData.getDoubleVector("vertical_coordinate_reference_pressure");
    }

    /// Units for a field, taken from the field metadata so written files carry the same units
    /// the rest of i-jedi believes the field is in. Empty when the name is not in the metadata.
    std::string unitsForField(const Geometry & geom, const std::string & longName)
    {
      const std::vector<std::string> & longNames = geom.getFieldMetadata().getLongNames();
      if (std::find(longNames.begin(), longNames.end(), longName) == longNames.end()) return "";
      return geom.getFieldMetadata().getFieldMetadata(longName).getVarUnits();
    }

    /// Define a CF coordinate variable with the attributes a reader expects to find.
    void defineCoord(int ncid, const std::string & name, int dimid, const std::string & units,
                     const std::string & standardName, const std::string & axis, int * varid)
    {
      ncCheck(nc_def_var(ncid, name.c_str(), NC_DOUBLE, 1, &dimid, varid),
              "defining coordinate " + name);
      nc_put_att_text(ncid, *varid, "units", units.size(), units.c_str());
      nc_put_att_text(ncid, *varid, "standard_name", standardName.size(), standardName.c_str());
      nc_put_att_text(ncid, *varid, "axis", axis.size(), axis.c_str());
    }

  }  // namespace

  // -------------------------------------------------------------------------------------------------

  IoStructuredNetcdf::IoStructuredNetcdf(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()), geometry_(geom), params_(params)
  {
    const std::string & precision = params_.precision.value();
    if (precision != "float" && precision != "double")
    {
      throw eckit::BadValue("IoStructuredNetcdf: 'precision' must be 'float' or 'double'",
                            Here());
    }
    const std::string & levelUnits = params_.levelUnits.value();
    if (levelUnits != "hPa" && levelUnits != "Pa")
    {
      throw eckit::BadValue("IoStructuredNetcdf: 'level units' must be 'hPa' or 'Pa'", Here());
    }
  }

  // -------------------------------------------------------------------------------------------------

  void IoStructuredNetcdf::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                                const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "read");
    oops::Log::trace() << classname() << " read starting" << std::endl;

    const std::string filepath = params_.filepath.value();
    const StructuredGridLayout layout = StructuredGridLayout::create(geometry_.functionSpace());
    const std::vector<double> geomLevelsPa = geometryLevelsPa(geometry_);

    // Every task opens the file for itself and reads only the block it owns. Reads are
    // independent, so this needs no parallel netCDF and no collectives.
    int ncid = 0;
    ncCheck(nc_open(filepath.c_str(), NC_NOWRITE, &ncid), "opening " + filepath);

    try
    {
      const FileMapping mapping =
          matchFile(ncid, layout, geomLevelsPa, params_.timeIndex.value(), filepath);

      int lonDimId = 0, latDimId = 0, levelDimId = -1, timeDimId = -1;
      ncCheck(nc_inq_dimid(ncid, mapping.lonDim.c_str(), &lonDimId), "longitude dimension");
      ncCheck(nc_inq_dimid(ncid, mapping.latDim.c_str(), &latDimId), "latitude dimension");
      if (!mapping.levelDim.empty()) nc_inq_dimid(ncid, mapping.levelDim.c_str(), &levelDimId);
      if (mapping.hasTime) nc_inq_dimid(ncid, mapping.timeDim.c_str(), &timeDimId);

      // First file row of this task's block, accounting for a file whose latitudes run the
      // other way round from the grid's.
      const int j0File = mapping.latFlipped ? (layout.ny - layout.j1) : layout.j0;
      const size_t boxNy = static_cast<size_t>(layout.boxNy());
      const size_t boxNx = static_cast<size_t>(layout.boxNx());

      for (atlas::Field & field : x)
      {
        const std::string jediName = field.name();
        const std::string fileName = fileionames.has(jediName)
                                     ? fileionames.getString(jediName) : jediName;
        const double scale = fileioscaling.has(jediName)
                             ? fileioscaling.getDouble(jediName) : 1.0;
        const int nLevels = field.shape(1);

        int varid = 0;
        ncCheck(nc_inq_varid(ncid, fileName.c_str(), &varid),
                "looking for variable '" + fileName + "' (JEDI name '" + jediName + "') in " +
                filepath);

        // Work out the hyperslab from the variable's own dimension list, so a file with or
        // without a time or level axis is handled without special-casing the caller.
        int nDims = 0;
        ncCheck(nc_inq_varndims(ncid, varid, &nDims), "rank of variable " + fileName);
        std::vector<int> dimIds(nDims);
        ncCheck(nc_inq_vardimid(ncid, varid, dimIds.data()), "dimensions of variable " + fileName);

        std::vector<size_t> start(nDims, 0);
        std::vector<size_t> count(nDims, 1);
        bool varHasLevel = false;
        size_t nLevelInVar = 1;

        for (int d = 0; d < nDims; ++d)
        {
          if (dimIds[d] == lonDimId)
          {
            start[d] = static_cast<size_t>(layout.i0);
            count[d] = boxNx;
            if (d != nDims - 1)
            {
              throw eckit::BadValue("IoStructuredNetcdf: variable '" + fileName + "' in " +
                                    filepath + " does not have longitude as its fastest-varying "
                                    "dimension", Here());
            }
          } else if (dimIds[d] == latDimId) {
            start[d] = static_cast<size_t>(j0File);
            count[d] = boxNy;
            if (d != nDims - 2)
            {
              throw eckit::BadValue("IoStructuredNetcdf: variable '" + fileName + "' in " +
                                    filepath + " does not have latitude immediately outside "
                                    "longitude", Here());
            }
          } else if (levelDimId >= 0 && dimIds[d] == levelDimId) {
            varHasLevel = true;
            nLevelInVar = mapping.nLevelFile;
            start[d] = 0;
            count[d] = nLevelInVar;
          } else if (timeDimId >= 0 && dimIds[d] == timeDimId) {
            start[d] = static_cast<size_t>(params_.timeIndex.value());
            count[d] = 1;
          } else {
            char dimName[NC_MAX_NAME + 1] = {0};
            nc_inq_dimname(ncid, dimIds[d], dimName);
            throw eckit::BadValue("IoStructuredNetcdf: variable '" + fileName + "' in " +
                                  filepath + " has unexpected dimension '" +
                                  std::string(dimName) + "'", Here());
          }
        }

        if (nLevels > 1)
        {
          if (geomLevelsPa.empty())
          {
            throw eckit::BadValue("IoStructuredNetcdf: reading the multi-level field '" +
                                  jediName + "' needs the geometry to carry a reference "
                                  "pressure column; set 'vertical levels' on the geometry",
                                  Here());
          }
          if (!varHasLevel)
          {
            throw eckit::BadValue("IoStructuredNetcdf: '" + jediName + "' needs " +
                                  std::to_string(nLevels) + " levels but '" + fileName +
                                  "' in " + filepath + " has no level dimension", Here());
          }
          if (static_cast<size_t>(nLevels) != mapping.levelIndex.size())
          {
            throw eckit::BadValue("IoStructuredNetcdf: '" + jediName + "' has " +
                                  std::to_string(nLevels) + " levels, which is not the "
                                  "geometry's level count; this backend handles full-level "
                                  "fields and single-level fields only", Here());
          }
        }

        if (boxNy == 0 || boxNx == 0) continue;  // task owns nothing; nothing to read

        // One call per field: the whole local block, all levels at once.
        std::vector<double> buffer(nLevelInVar * boxNy * boxNx);
        ncCheck(nc_get_vara_double(ncid, varid, start.data(), count.data(), buffer.data()),
                "reading variable " + fileName + " from " + filepath);

        // Scatter over the points this task actually owns. With more than a handful of tasks
        // atlas partitions in two dimensions, so the owned longitudes vary with latitude and
        // are a strict subset of the bounding box that was read; index(i, j) is meaningful
        // only inside them.
        auto view = atlas::array::make_view<double, 2>(field);
        for (int j = layout.j0; j < layout.j1; ++j)
        {
          const int fileRow = mapping.latFlipped ? (layout.ny - 1 - j) : j;
          const size_t r = static_cast<size_t>(fileRow - j0File);
          const int iOwnedBegin = static_cast<int>(layout.fs.i_begin(j));
          const int iOwnedEnd = static_cast<int>(layout.fs.i_end(j));
          for (int i = iOwnedBegin; i < iOwnedEnd; ++i)
          {
            const size_t c = static_cast<size_t>(i - layout.i0);
            const atlas::idx_t node = layout.index(i, j);
            for (int k = 0; k < nLevels; ++k)
            {
              const size_t kf = varHasLevel && nLevels > 1
                                ? static_cast<size_t>(mapping.levelIndex[k]) : 0;
              view(node, k) = buffer[(kf * boxNy + r) * boxNx + c] * scale;
            }
          }
        }
      }
    }
    catch (...)
    {
      nc_close(ncid);
      throw;
    }
    nc_close(ncid);

    // Owned points are filled above; fill the halos so downstream interpolation sees a
    // consistent field.
    geometry_.functionSpace().haloExchange(x);

    oops::Log::trace() << classname() << " read done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------

  void IoStructuredNetcdf::write(const atlas::FieldSet &x,
                                 const eckit::LocalConfiguration &fileionames,
                                 const eckit::LocalConfiguration &fileioscaling) const
  {
    util::Timer timer(classname(), "write");
    oops::Log::trace() << classname() << " write starting" << std::endl;

    const std::string filepath = params_.filepath.value();
    const StructuredGridLayout layout = StructuredGridLayout::create(geometry_.functionSpace());
    const std::vector<double> geomLevelsPa = geometryLevelsPa(geometry_);
    const eckit::mpi::Comm & comm = geometry_.getComm();
    const bool isRoot = (comm.rank() == 0);

    // Gather to the root task and write from there. At the resolutions the AI models run at a
    // global field is a few MB, so a single writer is not a bottleneck and avoids depending on
    // a parallel-enabled netCDF build.
    std::vector<atlas::Field> globalFields;
    std::vector<std::string> jediNames;
    for (const atlas::Field & field : x)
    {
      atlas::Field global = layout.fs.createField<double>(
          atlas::option::name(field.name()) | atlas::option::levels(field.shape(1)) |
          atlas::option::global());
      layout.fs.gather(field, global);
      globalFields.push_back(global);
      jediNames.push_back(field.name());
    }

    if (!isRoot)
    {
      oops::Log::trace() << classname() << " write done (non-root)" << std::endl;
      return;
    }

    const bool hasTime = params_.datetime.value() != boost::none;
    const bool levelsInHectoPascal = (params_.levelUnits.value() == "hPa");
    const bool asFloat = (params_.precision.value() == "float");

    int ncid = 0;
    ncCheck(nc_create(filepath.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid),
            "creating " + filepath);

    // Any field with more than one level needs the level axis.
    bool needLevelDim = false;
    for (const atlas::Field & field : globalFields)
    {
      if (field.shape(1) > 1) needLevelDim = true;
    }

    int lonDimId = 0, latDimId = 0, levelDimId = -1, timeDimId = -1;
    ncCheck(nc_def_dim(ncid, "longitude", layout.nx, &lonDimId), "defining longitude");
    ncCheck(nc_def_dim(ncid, "latitude", layout.ny, &latDimId), "defining latitude");
    if (needLevelDim)
    {
      ncCheck(nc_def_dim(ncid, "level", geomLevelsPa.size(), &levelDimId), "defining level");
    }
    if (hasTime)
    {
      ncCheck(nc_def_dim(ncid, "time", 1, &timeDimId), "defining time");
    }

    int lonVarId = 0, latVarId = 0, levelVarId = -1, timeVarId = -1;
    defineCoord(ncid, "longitude", lonDimId, "degrees_east", "longitude", "X", &lonVarId);
    defineCoord(ncid, "latitude", latDimId, "degrees_north", "latitude", "Y", &latVarId);
    if (needLevelDim)
    {
      defineCoord(ncid, "level", levelDimId, levelsInHectoPascal ? "hPa" : "Pa",
                  "air_pressure", "Z", &levelVarId);
    }
    if (hasTime)
    {
      defineCoord(ncid, "time", timeDimId, "hours since 1970-01-01T00:00:00", "time", "T",
                  &timeVarId);
    }

    // Data variables, in the canonical (time, level, latitude, longitude) order.
    std::vector<int> varIds(globalFields.size(), 0);
    for (size_t f = 0; f < globalFields.size(); ++f)
    {
      const std::string & jediName = jediNames[f];
      const std::string fileName = fileionames.has(jediName)
                                   ? fileionames.getString(jediName) : jediName;
      const bool is3D = globalFields[f].shape(1) > 1;

      std::vector<int> dims;
      if (hasTime) dims.push_back(timeDimId);
      if (is3D) dims.push_back(levelDimId);
      dims.push_back(latDimId);
      dims.push_back(lonDimId);

      ncCheck(nc_def_var(ncid, fileName.c_str(), asFloat ? NC_FLOAT : NC_DOUBLE,
                         static_cast<int>(dims.size()), dims.data(), &varIds[f]),
              "defining variable " + fileName);

      const std::string units = unitsForField(geometry_, jediName);
      if (!units.empty())
      {
        nc_put_att_text(ncid, varIds[f], "units", units.size(), units.c_str());
      }
    }

    ncCheck(nc_enddef(ncid), "leaving define mode for " + filepath);

    // Latitude order on disk. ARCO-ERA5 and WeatherBench-2 run south to north, while the
    // Gaussian grids atlas builds run north to south, so this is usually a flip.
    const bool flip = (params_.latStoN.value() == (layout.lats.front() > layout.lats.back()));

    std::vector<double> fileLats(layout.ny);
    for (int j = 0; j < layout.ny; ++j)
    {
      fileLats[j] = layout.lats[flip ? (layout.ny - 1 - j) : j];
    }
    ncCheck(nc_put_var_double(ncid, lonVarId, layout.lons.data()), "writing longitude");
    ncCheck(nc_put_var_double(ncid, latVarId, fileLats.data()), "writing latitude");
    if (needLevelDim)
    {
      std::vector<double> levelsOut(geomLevelsPa.size());
      for (size_t k = 0; k < geomLevelsPa.size(); ++k)
      {
        levelsOut[k] = levelsInHectoPascal ? geomLevelsPa[k] / 100.0 : geomLevelsPa[k];
      }
      ncCheck(nc_put_var_double(ncid, levelVarId, levelsOut.data()), "writing level");
    }
    if (hasTime)
    {
      const util::DateTime epoch("1970-01-01T00:00:00Z");
      const util::DateTime valid(params_.datetime.value().value());
      const double hours = static_cast<double>((valid - epoch).toSeconds()) / 3600.0;
      ncCheck(nc_put_var_double(ncid, timeVarId, &hours), "writing time");
    }

    const size_t nxs = static_cast<size_t>(layout.nx);
    const size_t nys = static_cast<size_t>(layout.ny);
    for (size_t f = 0; f < globalFields.size(); ++f)
    {
      const double scale = fileioscaling.has(jediNames[f])
                           ? fileioscaling.getDouble(jediNames[f]) : 1.0;
      const int nLevels = globalFields[f].shape(1);
      const auto view = atlas::array::make_view<double, 2>(globalFields[f]);

      std::vector<double> buffer(static_cast<size_t>(nLevels) * nys * nxs);
      for (int k = 0; k < nLevels; ++k)
      {
        for (int j = 0; j < layout.ny; ++j)
        {
          const int gridRow = flip ? (layout.ny - 1 - j) : j;
          for (int i = 0; i < layout.nx; ++i)
          {
            // The global field follows the grid's own row-major (j, i) ordering.
            const size_t src = static_cast<size_t>(gridRow) * nxs + static_cast<size_t>(i);
            buffer[(static_cast<size_t>(k) * nys + static_cast<size_t>(j)) * nxs +
                   static_cast<size_t>(i)] = view(src, k) / scale;
          }
        }
      }
      ncCheck(nc_put_var_double(ncid, varIds[f], buffer.data()),
              "writing variable " + jediNames[f]);
    }

    ncCheck(nc_close(ncid), "closing " + filepath);
    oops::Log::info() << "IoStructuredNetcdf: wrote " << filepath << std::endl;
    oops::Log::trace() << classname() << " write done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------

  void IoStructuredNetcdf::print(std::ostream &os) const
  {
    os << classname() << ": structured lat-lon netCDF (ERA5 / WeatherBench-2 layout), file "
       << params_.filepath.value();
  }

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

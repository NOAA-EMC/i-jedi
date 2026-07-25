/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <netcdf.h>

#include <map>
#include <algorithm>
#include <cmath>
#include <ostream>
#include <string>
#include <vector>


#include "atlas/grid.h"
#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/stringFunctions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
//#include "oops/base/GeometryData.h"

#include "ijedi/FieldMetadata/FieldsMetadata.h"
#include "ijedi/Increment/Increment.h"
#include "ijedi/State/State.h"

//#include "fv3jedi/Utilities/fv3jedi_vertical_remap.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/landVector/IOLandVector.h"



namespace ijedi {
// -------------------------------------------------------------------------------------------------
static IoMaker<IOLandVector> makerIOLandVector_("land vector");
// -------------------------------------------------------------------------------------------------
static inline void nc_rc(const int return_code, const std::string & operation) {
  if (return_code != NC_NOERR) {
    ABORT("IOLandVector netCDF operation \'" + operation + "\' failed with error: "
          + nc_strerror(return_code));
  }
}

IOLandVector::IOLandVector(const Geometry & geom, const Parameters_ & params) 
  : IoBase(geom, params.toConfiguration()), geom_(geom), params_(params) {
    
    // 1. All ranks participate natively in the parallel layout built during initialization
    // 2. No `writeFunctionSpace_` or `interpolator_` members needed!
    oops::Log::trace() << classname() << " Constructor configured for land vector." << std::endl;
}
// -------------------------------------------------------------------------------------------------
IOLandVector::~IOLandVector() {
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}

void IOLandVector::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                            const eckit::LocalConfiguration &fileioscaling) const 
{
  util::Timer timer(classname(), "read");
  oops::Log::trace() << classname() << " read started" << std::endl;

  // 1. Get your parallel FunctionSpace handle
  auto readFunctionSpace_ = atlas::functionspace::PointCloud(geom_.functionSpace());

  // 2. Create a temporary serial distribution where ALL points live on Rank 0
  //std::vector<int> zeros(geom_.globalNodeCount(), 0);
  //atlas::grid::Distribution serialDist(geom_.getComm().size(), geom_.globalNodeCount(), zeros.data());  
  const atlas::Grid & grid = readFunctionSpace_.grid();
  std::vector<int> zeros(grid.size(), 0);
  atlas::grid::Distribution serialDist(geom_.getComm().size(), grid.size(), zeros.data());

  eckit::LocalConfiguration atlas_conf;
  atlas_conf.set("mpi_comm", geom_.getComm().name());
  atlas::functionspace::PointCloud serialFunctionSpace(grid, serialDist, atlas_conf);

  // 3. Allocate fieldsSerial using the serial function space factory
  // This guarantees all ranks have perfectly matching metadata, ranks, and levels
  atlas::FieldSet fieldsSerial;
  for (const auto & field : x) {
    // Pass the matching metadata/options from the target field
    atlas::Field fSerial = serialFunctionSpace.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(field.levels()));
    fieldsSerial.add(fSerial);
  }

  int layer_index = 0;

  // 4. Loop through and read filenames on Rank 0
  const auto & filenamesOpt = params_.filenames.value();
  if (filenamesOpt != boost::none) {
    for (const auto & filename : filenamesOpt.value()) {
      if (geom_.getComm().rank() == 0) {
        util::DateTime dummyTime; 
        this->readVectorFields(params_.datapath.value() + "/" + filename, 
                                   fieldsSerial, dummyTime, layer_index, fileionames, fileioscaling);
      }
    }
  } else {
    ABORT("IOLandVector::read: 'filenames' parameter is missing.");
  }

  // 5. Scatter the data globally
  // Since both function spaces share the exact same grid, this will now succeed seamlessly
  readFunctionSpace_.scatter(fieldsSerial, x);
  
  oops::Log::trace() << classname() << " read done" << std::endl;
}


// -------------------------------------------------------------------------------------------------

void IOLandVector::readVectorFields(const std::string pathFile,
                                            atlas::FieldSet & fields,
                                            const util::DateTime & time,
                                            int layer_index = -1,
                                            const eckit::LocalConfiguration & ioNames,
                                            const eckit::LocalConfiguration & ioScaling) const {
  // NetCDF IDs
  int fileId;

  // Open a file to read fields from
  // -------------------------------
  oops::Log::info() << "Reading file " << pathFile << std::endl;
  nc_rc(nc_open(pathFile.c_str(), NC_NOWRITE, &fileId), "nc_open " + pathFile);

  int num_points = geom_.globalNodeCount();

  // Get file number of dimensions + their IDs
  // -----------------------------------------
  int ndims;
  // soil fields in form: soil_liquid_vol(time, soil_levels, location) ;
  nc_rc(nc_inq_ndims(fileId, &ndims), "nc_inq_ndims");

  std::vector<int> dimids(ndims);
  nc_rc(nc_inq_dimids(fileId, &ndims, dimids.data(), 0), "nc_inq_dimids");

  
  // Ensure that the lat and lon dimensions are found and have the correct lengths
  size_t dimSize;
  bool hasLoc = false;
  bool hasLayer = false;
  bool hasTim = false;
  int locId;
  int layerId;
  int timId;
  for (int i = 0; i < ndims; ++i) {
    // Get the name and size of the dimension
    char dimName[NC_MAX_NAME + 1];
    nc_rc(nc_inq_dim(fileId, dimids[i], dimName, &dimSize), "nc_inq_dim");

    if (std::string(dimName) == params_.locatoinName.value().c_str()) {
      hasLoc = true;
      locId = dimids[i];
      ASSERT(dimSize == num_points);
    } else if (std::string(dimName) == params_.layerName.value().c_str()) {
      hasLayer = true;
      layerId = dimids[i];
      ASSERT(dimSize > layer_index);  //TODO: Check if this is the correct assertion for layer_index
    } else if (std::string(dimName) == params_.timName.value().c_str()) {
      hasTim = true;
      timId = dimids[i];
      ASSERT(dimSize == 1); // Only one time step is expected for this read
    }
  }

  // Ensure required dimensions were found
  ASSERT(hasLoc);
  ASSERT(hasTim);
  ASSERT(hasLayer || layer_index == -1); // no layer dimension
  
  // Read the fields from the file
  // -------------------------------
  for (auto & field : fields) {
    // Get IO name for this field
    std::string fieldName = field.name();
    if (ioNames.has(fieldName)) {
      fieldName = ioNames.getString(field.name());
    }
    oops::Log::info() << "Field " << fieldName << std::endl;
    // Get the variable ID for this field
    int varId;
    int status = nc_inq_varid(fileId, fieldName.c_str(), &varId);
    if (status == NC_ENOTVAR) {
      // Variable is not in this file; skip it silently (or log an info message)
      oops::Log::info() << "Field " << fieldName << " not found in this file, skipping..." << std::endl;
      continue; 
    } else {
      // Check for any other unexpected NetCDF errors
      nc_rc(status, "nc_inq_varid " + fieldName);
    }
    // nc_rc(nc_inq_varid(fileId, fieldName.c_str(), &varId), "nc_inq_varid " + fieldName);


    // Get number of dimensions + their IDs
    int vardimids[NC_MAX_VAR_DIMS];
    nc_rc(nc_inq_var(fileId, varId,
                     nullptr,   // var name (unused)
                     nullptr,   // type (unused)
                     &ndims,
                     vardimids,
                     nullptr),  // attributes (unused)
          "nc_inq_var");

    ASSERT(ndims == 2 || ndims == 3);

    // Ensure that the dimensions are in the expected order
    if ( ndims == 2 ) {
      ASSERT(vardimids[0] == timId && vardimids[1] == locId);
      size_t start[2] = {0, 0};
      size_t count[2] = {1, num_points};

    } else if ( ndims == 3 ) {
      ASSERT(vardimids[0] == timId && vardimids[1] == layerId && vardimids[2] == locId);
      size_t start[3] = {0, layer_index, 0};
      size_t count[3] = {1, 1, num_points};
    }

    // Read the variable data
    std::vector<double> values(field.size());

    nc_rc(nc_get_vara_double(fileId, varId, start, count, values.data()), "nc_get_var_double " + fieldName);

    // Create field and unpack data into it
    /*if (field.rank() == 2) {
      // Standard multi-level or Rank-2 surface field [Points, layers]
      auto fieldView = atlas::array::make_view<double, 2>(field);

      for (size_t k = 0; k < numLayers; ++k) {
          for (size_t i = 0; i < num_points; ++i) {
            fieldView(i, k) = values[ k*num_points + i ];
          }
        }
      }
    } else */
    if (field.rank() == 1) {
      // Pure Rank-1 surface field [Points]
      auto fieldView = atlas::array::make_view<double, 1>(field);
      ASSERT(numLayers == 1);
      for (size_t j = 0; j < num_points; ++j) {
          fieldView(j) = values[j];
      }
    } else {
      ABORT("IOLandVector::readVectorFields - Unsupported field rank: " + std::to_string(field.rank()));
    }
    oops::Log::info() << "Done reading Field " << fieldName << std::endl;
  }
  // Close file
  nc_rc(nc_close(fileId), "nc_close");
}

void IOLandVector::write(const atlas::FieldSet & fieldsVector,
                             const eckit::LocalConfiguration & fileionames,
                             const eckit::LocalConfiguration & fileioscaling) const 
{

  util::Timer timer(classname(), "write");
  oops::Log::trace() << classname() << " write started" << std::endl;

  // 1. Get your parallel FunctionSpace handle
  auto readFunctionSpace_ = atlas::functionspace::PointCloud(geom_.functionSpace());

  // 2. Create a temporary serial distribution where ALL points live on Rank 0
  //std::vector<int> zeros(geom_.globalNodeCount(), 0);
  //atlas::grid::Distribution serialDist(geom_.getComm().size(), geom_.globalNodeCount(), zeros.data());  
  const atlas::Grid & grid = readFunctionSpace_.grid();
  std::vector<int> zeros(grid.size(), 0);
  atlas::grid::Distribution serialDist(geom_.getComm().size(), grid.size(), zeros.data());

  eckit::LocalConfiguration atlas_conf;
  atlas_conf.set("mpi_comm", geom_.getComm().name());
  atlas::functionspace::PointCloud serialFunctionSpace(grid, serialDist, atlas_conf);

  // 3. Allocate fieldsSerial using the serial function space factory
  // This guarantees all ranks have perfectly matching metadata, ranks, and levels
  atlas::FieldSet fieldsSerial;
  for (const auto & field : x) {
    // Pass the matching metadata/options from the target field
    atlas::Field fSerial = serialFunctionSpace.createField<double>(
        atlas::option::name(field.name()) | atlas::option::levels(field.levels()));
    fieldsSerial.add(fSerial);
  }

  // Gather distributed data smoothly from all parallel ranks back into Rank 0
  readFunctionSpace_.gather(fieldsVector, fieldsSerial);

  // Resolve the valid time dynamically
  util::DateTime validTime;
  if (fieldsVector.metadata().has("time")) {
    validTime = util::DateTime(fieldsVector.metadata().get<std::string>("time"));
  } else if (fieldsVector.metadata().has("timestamp")) {
    validTime = util::DateTime(fieldsVector.metadata().get<std::string>("timestamp"));
  } else if (params_.dateTime.value() != boost::none) {
    validTime = util::DateTime(params_.dateTime.value().value());
  } else {
    oops::Log::warning() << "Using default DateTime placeholder for file writing." << std::endl;
    validTime = util::DateTime();
  }

  // Write to disk exclusively on rank 0
  if (geom_.getComm().rank() == 0) {
    //const util::DateTime dateTime(datTimeString);
    this->writeVectorFields(fieldsSerial, validTime, fileionames, fileioscaling);
  }
  oops::Log::trace() << classname() << " write done" << std::endl;

}

void IOLandVector::writeVectorFields(const atlas::FieldSet & fields,
                                             const util::DateTime & time,
                                             const eckit::LocalConfiguration & ioNames,
                                             const eckit::LocalConfiguration & ioScaling) const {
  
  // NetCDF IDs
  // ----------
  int fileId, fIv, locId, layerId, timId;
  std::map<std::string, int> fieldIvs;
  int nTim = 1;

  // Get the name of the file and adjust with datetime
  // -------------------------------------------------
  std::string pathFile = params_.filename.value();

  // Replace member number (ensemble applciaitons)
  util::stringfunctions::swapNameMember(params_.toConfiguration(), pathFile);

  // Create a file to write fields into
  // ----------------------------------
  nc_rc(nc_create(pathFile.c_str(), NC_CLOBBER | NC_NETCDF4, &fileId), "nc_create" + pathFile);
  oops::Log::warning() << "nc created" << std::endl;

  int num_locations = geom_.globalNodeCount();

  // Set float precision for fields
  // ------------------------------
  const int floatPrecision = params_.floatPrecision.value();
  const int ncPrec = (floatPrecision == 4) ? NC_FLOAT : NC_DOUBLE;

  //Get time in seconds since epoch
  const util::DateTime epoch("1970-01-01T00:00:00Z");
  const util::Duration duration = time - epoch;
  int64_t seconds_since_epoch = duration.toSeconds();


  nc_rc(nc_def_dim(fileId, params_.locationName.value().c_str(), num_locations, &locId), "nc_def_dim (location)");
  //nc_rc(nc_def_dim(fileId, params_.layerName.value().c_str(), num_layers, &layerId), "nc_def_dim (layer)");
  nc_rc(nc_def_dim(fileId, params_.timeName.value().c_str(), nTim, &timId), "nc_def_dim (time)");

  // Write the dimension variables: only time for now
  nc_rc(nc_def_var(fileId, params_.timeName.value().c_str(), NC_INT, 1, &timId, &fIv),
        "nc_def_var (tim)");
  nc_rc(nc_put_att_text(fileId, fIv, "long name", strlen("time"), "time"),
        "nc_put_att_text (long name time)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("seconds since 1970-01-01 00:00:00"), "seconds since 1970-01-01 00:00:00"),
        "nc_put_att_text (units time)");
  fieldIvs[params_.timeName.value()] = fIv;

  // Define some categories of dimension IDs for fields
  // --------------------------------------------------
  std::vector<int> fieldDims = {timId, locId};  // only one layer written out

  // Define all the fields that will be written
  // ------------------------------------------
  for (auto& field : fields) {
    
    ASSERT(field.shape(1) == num_locations);  // Ensure the field has the expected number of locations

    // Get dimensions for this field 
    const auto &dims = field.shape();

    // Look for fieldname in the iofile configuration and use the value if key found
    const std::string fieldLong = field.name();
    const char * fieldLongC = fieldLong.c_str();
    
    std::string fieldName = fieldLong;
    if (ioNames.has(fieldName)) {
      fieldName = ioNames.getString(fieldLong);
    }

    // Define the field in the file
    nc_rc(nc_def_var(fileId, fieldName.c_str(), ncPrec, dims.size(), dims.data(), &fIv), "nc_def_var " + fieldName);

    // Fallback defaults if metadata keys are missing
    std::string unitsStr = "unknown";
    std::string longNameStr = fieldLong;

    // Extract values dynamically if Atlas has them populated
    if (field.metadata().has("units")) {
      unitsStr = field.metadata().get<std::string>("units");
    }
    if (field.metadata().has("long_name")) {
      longNameStr = field.metadata().get<std::string>("long_name");
    }

    const char * units = unitsStr.c_str();
    const char * longNameC = longNameStr.c_str();

    // Write to NetCDF
    nc_rc(nc_put_att_text(fileId, fIv, "units", strlen(units), units), "nc_put_att_text " + fieldName + " units");
    nc_rc(nc_put_att_text(fileId, fIv, "long_name", strlen(longNameC), longNameC), 
          "nc_put_att_text " + fieldName + " long_name");

    // Insert field into the fieldIvs map
    fieldIvs[field.name()] = fIv;
  }
  // End definition mode
  // -------------------
  nc_rc(nc_enddef(fileId), "nc_enddef");

  // Write coordinate data 
  nc_rc(nc_put_var_int(fileId, fieldIvs[params_.timeName.value()], seconds_since_epoch), "nc_put_var_int (time)");

  // Write the fields into the file
  // ------------------------------
  for (auto& field : fields) {
 
    // Create view of the field
    const auto fieldView = atlas::array::make_view<double, 1>(field);

    // Vector to hold the packed field
    std::vector<double> values(num_locations);

    // Loop over dimensions and pack the field
    for (size_t k = 0; k < num_locations; ++k) {
          values[k] = fieldView(k);
    }

    // Write the field to the file
    nc_rc(nc_put_var_double(fileId, fieldIvs[field.name()], values.data()), "nc_put_var_double " + field.name());
  }

  // Close netCDF file
  // -----------------
  nc_rc(nc_close(fileId), "nc_close");
}

// -------------------------------------------------------------------------------------------------

void IOLandVector::print(std::ostream & os) const {
  os << classname() << " IO for land vector using Atlas PointCloud FunctionSpace";
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jed

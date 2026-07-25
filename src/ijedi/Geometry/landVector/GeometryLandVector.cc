/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

 // A geometry class that only needs lat/lon/elevation information--for land (snow/soil) DA over independent grids.
 // Written with the intent of being a minimal geometry class that can be used for land DA over independent grids, without requiring a full model geometry.
 // Written with the help of Gemeni, fully reviewed and tested by land DA team.
  

#include "ijedi/Geometry/landVector/GeometryLandVector.h"

#include "atlas/array.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/util/Config.h"
#include "oops/util/Logger.h"

namespace ijedi {

// -----------------------------------------------------------------------------
GeometryLandVector::GeometryLandVector(const eckit::Configuration &config,
                                       const eckit::mpi::Comm &comm,
                                       eckit::LocalConfiguration & /*geomVariables*/,
                                       atlas::FunctionSpace &functionSpace,
                                       atlas::FieldSet &geomFields,
                                       bool &levelsAreTopDown, int &numLevels
    : comm_(comm) {
  oops::Log::trace() << "GeometryLandVector::GeometryLandVector starting..." << std::endl;
 
  numLevels_ = config.getInt("nlevels", 1);

  // 1. Read coordinates (lat, lon, elevation) from configuration or input file
  std::string latlonFile = config.getString("latlon_file", "");
  std::string latName = config.getString("lat_name", "latitude");
  std::string lonName = config.getString("lon_name", "longitude");
  std::string elevName = config.getString("elev_name", "elevation");

  if(!latlonFile.empty() && config.has("lat_name") && config.has("lon_name") && config.has("elev_name")) {
    oops::Log::info() << "Reading lat/lon/elevation from file: " << latlonFile << std::endl;
    readLatLonElevFromFile(latlonFile, latName, lonName, elevName, lats, lons, elevations);
  } else if (config.has("lats") && config.has("lons") && config.has("elevations")) {
    oops::Log::info() << "Reading lat/lon/elevation from configuration." << std::endl;
    std::vector<double> lats = config.getDoubleVector("lats");
    std::vector<double> lons = config.getDoubleVector("lons");
    std::vector<double> elevations = config.getDoubleVector("elevations");
  }
  else {
    throw eckit::BadParameter("GeometryLandVector: Must specify either latlon_file or lats/lons/elevations in configuration.", Here());
  }

  ASSERT(lats.size() == lons.size());
  ASSERT(lats.size() == elevations.size());

  size_t numPoints = lats.size();

  // Global point cloud geometry
  std::vector<atlas::PointXY> global_pts(numPoints);
  for (size_t i = 0; i < numPoints; ++i) {
    global_pts[i] = atlas::PointXY(lons[i], lats[i]);
  }

  // global functionspace, for reading fields serially from file and then distributing to MPI tasks
  global_functionSpace_ = atlas::functionspace::PointCloud(global_pts);

  /*global_lonlatField_ = global_functionSpace_.createField<double>(
      atlas::option::name("lonlat") | atlas::option::variables(2));
  
  global_elevationField_ = global_functionSpace_.createField<double>(
      atlas::option::name("height") | atlas::option::variables(1));

  // Populate field values
  auto global_lonlatView = atlas::array::make_view<double, 2>(global_lonlatField_);
  auto global_elevView   = atlas::array::make_view<double, 1>(global_elevationField_);

  for (size_t i = 0; i < numPoints; ++i) {
    global_lonlatView(i, 0) = lons[i];
    global_lonlatView(i, 1) = lats[i];
    global_elevView(i)      = elevations[i];
  }*/

  atlas::UnstructuredGrid grid(global_pts);

  // Create partitioner (EqualRegions splits points evenly across MPI tasks)
  atlas::grid::Partitioner partitioner("equal_regions", comm_.size());
  atlas::grid::Distribution distribution = partitioner.partition(grid);

  // Extract points owned by rank
  int myRank = comm_.rank();
  std::vector<atlas::PointXY> local_pts;
  std::vector<double> local_elevations;

  for (size_t i = 0; i < numPoints; ++i) {
    if (distribution.partition(i) == myRank) {
      local_pts.push_back(global_pts[i]);
      local_elevations.push_back(elevations[i]);
    }
  }

  // Initialize PointCloud functionspace across local MPI tasks
  functionSpace_ = atlas::functionspace::PointCloud(local_pts);

  // Register standard JEDI coordinate fields (lonlat and vertical height/elevation)
  lonlatField_ = functionSpace_.createField<double>(
      atlas::option::name("lonlat") | atlas::option::variables(2));
  
  elevationField_ = functionSpace_.createField<double>(
      atlas::option::name("height") | atlas::option::variables(1));

  // Populate field values
  auto lonlatView = atlas::array::make_view<double, 2>(lonlatField_);
  auto elevView   = atlas::array::make_view<double, 1>(elevationField_);

  for (size_t i = 0; i < local_points.size(); ++i) {
    lonlatView(i, 0) = local_points[i].x();   //lons;
    lonlatView(i, 1) = local_points[i].y();   //lats;
    elevView(i)      = local_elevations[i];
  }

  // Add fields to geometry fieldset so JEDI ObsOperators / Increment objects can query them
  fields_.add(lonlatField_);
  fields_.add(elevationField_);

  oops::Log::trace() << "GeometryLandVector::GeometryLandVector completed with "
                     << numPoints << " global and " << local_points.size() << " local points." << std::endl;

}

// -----------------------------------------------------------------------------
GeometryLandVector::GeometryLandVector(const GeometryLandVector & other)
    : comm_(other.comm_),
      functionSpace_(other.functionSpace_),
      fields_(other.fields_),
      lonlatField_(other.lonlatField_),
      elevationField_(other.elevationField_) {}

// -----------------------------------------------------------------------------
GeometryLandVector::~GeometryLandVector() {}

// -----------------------------------------------------------------------------
size_t GeometryLandVector::globalNodeCount() const {
  return functionSpace_.globalSize();
}

// -----------------------------------------------------------------------------
size_t GeometryLandVector::localNodeCount() const {
  return functionSpace_.size();
}

// -----------------------------------------------------------------------------
void GeometryLandVector::print(std::ostream & os) const {
  os << "GeometryLandVector: [ Count = " << globalNodeCount() << " points ]";
}

// -----------------------------------------------------------------------------------------------

std::vector<double> GeometryLandVector::verticalCoord(std::string & /*vcUnits*/) const {
  std::stringstream errorMsg;
  errorMsg << "GeometryLandVector::verticalCoord is not implemented" << std::endl;
  throw eckit::NotImplemented(errorMsg.str(), Here());
}

void GeometryLandVector::readLatLonElevFromFile(const std::string pathFile, const std::string latName, const std::string lonName, const std::string elevName,
    std::vector<double> & lats, std::vector<double> & lons, std::vector<double> & elevations)
{
    // NetCDF IDs
    int fileId;

    oops::Log::info() << "Reading file " << pathFile << std::endl;
    nc_rc(nc_open(pathFile.c_str(), NC_NOWRITE, &fileId), "nc_open " + pathFile);

    // Get the variable ID for this field
    int varId;
    nc_rc(nc_inq_varid(fileId, latName.c_str(), &varId), "nc_inq_varid " + latName);

    // Get number of dimensions + their IDs
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    nc_rc(nc_inq_var(fileId, varId,
                     nullptr,   // var name (unused)
                     nullptr,   // type (unused)
                     &ndims,
                     dimids,
                     nullptr),  // attributes (unused)
          "nc_inq_var");

    // Ensure that we are working with a vector (1D variables)
    ASSERT(ndims == 1);
    // Read the variable data
    lats.resize(dimids[0]);
    lons.resize(dimids[0]);
    elevations.resize(dimids[0]);
    
    nc_rc(nc_get_var_double(fileId, varId, lats.data()), "nc_get_var_double " + latName);

    nc_rc(nc_inq_varid(fileId, lonName.c_str(), &varId), "nc_inq_varid " + lonName);
    nc_rc(nc_get_var_double(fileId, varId, lons.data()), "nc_get_var_double " + lonName);

    nc_rc(nc_inq_varid(fileId, elevName.c_str(), &varId), "nc_inq_varid " + elevName);
    nc_rc(nc_get_var_double(fileId, varId, elevations.data()), "nc_get_var_double " + elevName);

    oops::Log::info() << "Done reading Field " << pathFile << std::endl;
  
    // Close file
    nc_rc(nc_close(fileId), "nc_close "+ pathFile);
}

// -----------------------------------------------------------------------------------------------

}  // namespace ijedi
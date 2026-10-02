#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

namespace ijedi {
namespace mom6 {

int computeExtent(int N, int ndivs, int pe);

void readHgrid(const std::string & path,
               int niGlobal,
               int njGlobal,
               int niEff,
               int njEff,
               int coarsenFactor,
               std::vector<double> * lon,
               std::vector<double> * lat,
               std::vector<double> * dxT,
               std::vector<double> * dyT,
               std::vector<double> * areaT,
               std::vector<double> * lonU,
               std::vector<double> * latU,
               std::vector<double> * lonV,
               std::vector<double> * latV);

void readTopog(const std::string & path,
               int niGlobal,
               int njGlobal,
               int niEff,
               int njEff,
               int coarsenFactor,
               double minimumDepth,
               std::vector<double> * depth,
               std::vector<double> * wet);

void readVerticalGeometry(const std::string & path,
                         int niGlobal,
                         int njGlobal,
                         int niEff,
                         int njEff,
                         int numLevels,
                         int coarsenFactor,
                         std::vector<double> * layerThickness,
                         std::vector<double> * layerCenterDepth);

// Read 1-D point data (latitude, longitude, vars) from a netCDF file and interpolate
// the variables onto dstFunctionSpace with 10-nearest-neighbours weighting, as soca does
// for its geometry fields (e.g. rossby_radius from rossrad.nc).
atlas::FieldSet readNcAndInterp(const std::string & filename,
                                const std::vector<std::string> & vars,
                                const atlas::FunctionSpace & dstFunctionSpace);

void rewriteLeveledNodeDataAsTime(const std::string & filename,
                                  const std::unordered_set<std::string> & timeFields);

}  // namespace mom6
}  // namespace ijedi

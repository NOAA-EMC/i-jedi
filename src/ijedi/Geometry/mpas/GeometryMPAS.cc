#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

//#include "atlas/array/ArrayView.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/option.h"

//#include "oops/base/GeometryData.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ijedi/Geometry/mpas/GeometryMPAS.h"
#include "ijedi/Geometry/mpas/GeometryMPAS.interface.h"
//#include "ijedi/Geometry/mpas/GeometryMPASParameters.h"

namespace ijedi
{

  atlas::functionspace::NodeColumns buildNodeColumns(void * fortranGeom,
                                                     const eckit::mpi::Comm & comm) {

  // -----------------------------------------------------------------------------------------------

  GeometryMPAS::GeometryMPAS(const eckit::Configuration &geomConfig,
                             const eckit::mpi::Comm &comm,
                             eckit::Configuration &geomVariables,
                             atlas::FunctionSpace &functionSpace,
                             atlas::FieldSet &fieldSet,
                             int &numberLevels) 
  { 
  oops::Log::trace() << "GeometryMPAS constructor starting" << std::endl;

  //ijedi_mpas_geom_setup_f90(fortranGeom_, config, &comm);
  ijedi_mpas_geom_setup_f90(fortranGeom_, geomConfig, &comm);

  // Build Atlas NodeColumns function space from MPAS grid topology
  functionSpace = buildNodeColumns(fortranGeom_, comm);

  // Get nVertLevels from Fortran for vert_coord and levelsPerVariable_
  int nVertLevels;
  ijedi_mpas_geom_get_vertical_resolution_f90(fortranGeom_, nVertLevels);

  // Fill geometry fields (owned, area, vert_coord)
  fields_ = atlas::FieldSet();
  fillGeometryFields(fortranGeom_, functionspace_, nVertLevels, fields_);

  // Set vertical metadata
  levelsAreTopDown_  = true;
  levelsPerVariable_ = createLevelsPerVariable(nVertLevels);

  // Build GeometryData
  geomData_.reset(new oops::GeometryData(functionspace_, fields_, levelsAreTopDown_, comm));

  oops::Log::trace() << "ijedi_mpas::Geometry::Geometry from config done" << std::endl;


  }

  // -----------------------------------------------------------------------------------------------

  void GeometryMPAS::print(std::ostream &os) const 
  {
    int nVertLevels;
    ijedi_mpas_geom_get_vertical_resolution_f90(fortranGeom_, nVertLevels);
    int nCellsGlobal;
    ijedi_mpas_geom_get_global_cell_count_f90(fortranGeom_, nCellsGlobal);
    os << "ijedi_mpas::Geometry, nCellsGlobal = " << nCellsGlobal << ", nVertLevels = " << nVertLevels
       << ", communicator = " << comm().name();
  }

  // -----------------------------------------------------------------------------------------------

  std::vector<double> GeometryMPAS::verticalCoord(std::string &vcUnits) const
  {
    // Not implemented, abort --- IGNORE ---
    std::stringstream errorMsg;
    errorMsg << "GeometryMPAS::verticalCoord is not implemented" << std::endl;
    ABORT(errorMsg.str());
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace ijedi

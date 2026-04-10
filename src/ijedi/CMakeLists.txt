include_directories(${CMAKE_CURRENT_SOURCE_DIR})

list(APPEND ijedi_src_files
  Traits.h
)

# ###############################################################################
# Project source files
# ###############################################################################
list(APPEND ijedi_src_files

  # Geometry
  Geometry/Geometry.cc
  Geometry/Geometry.h
  Geometry/base/GeometryBase.cc
  Geometry/base/GeometryBase.h
  Geometry/mpas/GeometryMPASLoader.cc
  Geometry/mom6/GeometryMOM6.cc
  Geometry/mom6/GeometryMOM6.h

  # Fields
  FieldMetadata/FieldsMetadata.cc
  FieldMetadata/FieldsMetadata.h
  FieldMetadata/FieldsMetadataDefault.h

  # State
  State/State.cc
  State/State.h

  # Increment
  Increment/Increment.cc
  Increment/Increment.h

  # IO
  Io/IoBase.cc
  Io/IoBase.h
  Io/mom6/IoMOM6.cc
  Io/mom6/IoMOM6.h
  Io/mom6/ReadMOM6Netcdf.cc
  Io/mom6/ReadMOM6Netcdf.h
  Io/mom6/WriteMOM6Netcdf.cc
  Io/mom6/WriteMOM6Netcdf.h

  # Variable change
  VariableChange/SeaSurfaceTemperature.h
  VariableChange/SeaSurfaceTemperature_A.cc
  VariableChange/VaderCookbook.cc
  VariableChange/VaderCookbook.h
  VariableChange/VariableChange.cc
  VariableChange/VariableChange.h

  # Linear variable change
  LinearVariableChange/LinearVariableChange.cc
  LinearVariableChange/LinearVariableChange.h

  # Model and LinearModel (persistence and identity)
  Model/Model.h
  LinearModel/LinearModel.h

  # Utilities
  Utilities/Constants.cc
  Utilities/Constants.h
  Utilities/ijedi_constants_mod.f90
  Utilities/ijedi_kinds_mod.f90
  Utilities/ijedi_netcdf_utils_mod.F90
  Utilities/PrintHelper.cc
  Utilities/PrintHelper.h
)

# ###############################################################################
# FV3 specific project source files and setup
# ###############################################################################

# fv3 pieces if FMS is found
if(FMS_FOUND)
  list(APPEND ijedi_fv3_cpp_src_files
    Geometry/fv3/GeometryFV3.cc
    Geometry/fv3/GeometryFV3.h
    Geometry/fv3/GeometryFV3Parameters.h
    Geometry/fv3/GeometryFV3.interface.h
    Io/fv3/IoFV3.cc
    Io/fv3/IoFV3.h
  )
  list(APPEND ijedi_fv3_fortran_src_files
    Geometry/fv3/ijedi_fv3_akbk_mod.f90
    Geometry/fv3/ijedi_fv3_arrays.F90
    Geometry/fv3/ijedi_fv3_control.F90
    Geometry/fv3/ijedi_fv3_eta.F90
    Geometry/fv3/ijedi_fv3_geom_interface_mod.f90
    Geometry/fv3/ijedi_fv3_geom_mod.f90
    Geometry/fv3/ijedi_fv3_grid_tools.F90
    Geometry/fv3/ijedi_fv3_grid_utils.F90
    Geometry/fv3/ijedi_fv3_mp_mod.F90
    Geometry/fv3/ijedi_fv3_namelist_mod.f90
    Geometry/fv3/ijedi_fv3_sorted_index.F90
  )

  list(APPEND ijedi_fv3_src_files
    ${ijedi_fv3_cpp_src_files}
    ${ijedi_fv3_fortran_src_files}
  )

  # fv3 ported files have specific required properties
  list(APPEND fv3_defs SPMD)
  list(APPEND fv3_defs MOIST_CAPPA)
  list(APPEND fv3_defs USE_COND)
  list(APPEND fv3_defs OVERLOAD_R4)
  list(APPEND fv3_defs OVERLOAD_R8)
  list(APPEND fv3_defs INTERNAL_FILE_NML)
  list(APPEND fv3_defs ENABLE_QUAD_PRECISION)
  set_source_files_properties(${ijedi_fv3_src_files} PROPERTIES COMPILE_DEFINITIONS "${fv3_defs}")

  # Combine with main source files
  list(APPEND ijedi_src_files ${ijedi_fv3_src_files})

  # Compile with double precision
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(fv3_flags "-r8")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(fv3_flags "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none")
  else()
    message(WARNING
      "Unknown Fortran compiler '${CMAKE_Fortran_COMPILER_ID}',\n"
      "  fv3 files not using special flags.")
    set(fv3_flags "")
  endif()

  set_source_files_properties(${ijedi_fv3_fortran_src_files} PROPERTIES COMPILE_FLAGS "${fv3_flags}")
endif()

# ###############################################################################
# MPAS specific project source files and setup
# ###############################################################################
if(BUILD_MPAS)

  list(APPEND ijedi_mpas_src_files
    Geometry/mpas/GeometryMPASPlugin.cc
    Geometry/mpas/GeometryMPAS.cc
    Geometry/mpas/GeometryMPAS.h
    Geometry/mpas/GeometryMPAS.interface.h
    Geometry/mpas/ijedi_mpas_geom_interface_mod.F90
    Geometry/mpas/ijedi_mpas_geom_mod.F90
  )

endif()

# ###############################################################################
# Model specific project source files
# ###############################################################################
ecbuild_add_library(TARGET ${PROJECT_NAME}
  SOURCES ${ijedi_src_files}
  INSTALL_HEADERS LISTED
  LINKER_LANGUAGE CXX
)

target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}>)

target_link_libraries(${PROJECT_NAME} PUBLIC NetCDF::NetCDF_C NetCDF::NetCDF_CXX)

# target_link_libraries( ${PROJECT_NAME} PUBLIC MPI::MPI_Fortran MPI::MPI_CXX )
target_link_libraries(${PROJECT_NAME} PUBLIC eckit)
target_link_libraries(${PROJECT_NAME} PUBLIC fckit)
target_link_libraries(${PROJECT_NAME} PUBLIC atlas atlas_f)
target_link_libraries(${PROJECT_NAME} PUBLIC oops)
target_link_libraries(${PROJECT_NAME} PUBLIC ${CMAKE_DL_LIBS})

# Mist
target_link_libraries(${PROJECT_NAME} PUBLIC mist_base)
target_link_libraries(${PROJECT_NAME} PUBLIC mist_utils)
# target_link_libraries( ${PROJECT_NAME} PUBLIC ufo )
if(FMS_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC fms)
endif()
if(BUILD_MPAS)
  ecbuild_add_library(TARGET ${PROJECT_NAME}_mpas
    SOURCES ${ijedi_mpas_src_files}
    INSTALL_HEADERS LISTED
    LINKER_LANGUAGE CXX
  )

  target_include_directories(${PROJECT_NAME}_mpas PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>)

  target_link_libraries(${PROJECT_NAME}_mpas PUBLIC ${PROJECT_NAME})
  target_link_libraries(${PROJECT_NAME}_mpas PUBLIC MPAS::core::atmosphere)
  target_compile_options(${PROJECT_NAME}_mpas PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

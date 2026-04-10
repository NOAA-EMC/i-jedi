# Build the test executable
# -------------------------
add_subdirectory(mains)

# Automatically discover and add tests organized by tier
# -------------------------------------------------------
# All tests are built; use ctest -R tier<N> to run specific tiers
# Note: CONFIGURE_DEPENDS tells CMake to recheck for new/removed test files
# `export IJEDI_MPI_ARGS="--map-by :OVERSUBSCRIBE"` to allow MPI oversubscription (useful for CI/local testing)
set(IJEDI_MPI_ARGS $ENV{IJEDI_MPI_ARGS})

if(IJEDI_MPI_ARGS)
    separate_arguments(IJEDI_MPI_ARGS)
endif()

# Link the directory 'data' to the test directory of the build
# ------------------------------------------------------------
file(CREATE_LINK
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    "${CMAKE_CURRENT_BINARY_DIR}/data"
    SYMBOLIC
)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data_generated)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data_generated/mom6)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data_generated/obs)

# Fetch MOM.res.nc from soca repo (sparse clone + LFS pull, skip if exists)
set(_MOM_RES_NC "${CMAKE_CURRENT_BINARY_DIR}/data_generated/mom6/MOM.res.nc")
if(NOT EXISTS "${_MOM_RES_NC}")
  message(STATUS "Fetching MOM.res.nc from JCSDA/soca ...")
  set(_SOCA_CLONE_DIR "${CMAKE_CURRENT_BINARY_DIR}/data_generated/mom6/_soca_sparse")
  file(REMOVE_RECURSE "${_SOCA_CLONE_DIR}")
  execute_process(
    COMMAND git clone --depth 1 --filter=blob:none --sparse
            https://github.com/JCSDA/soca.git "${_SOCA_CLONE_DIR}"
    RESULT_VARIABLE _clone_rc
    OUTPUT_QUIET)
  if(NOT _clone_rc EQUAL 0)
    message(FATAL_ERROR "Failed to sparse-clone JCSDA/soca")
  endif()
  execute_process(
    COMMAND git -C "${_SOCA_CLONE_DIR}" sparse-checkout set --no-cone /test/Data/72x35x25/restarts/MOM.res.nc
    RESULT_VARIABLE _sparse_rc
    OUTPUT_QUIET)
  if(NOT _sparse_rc EQUAL 0)
    message(FATAL_ERROR "Failed to set sparse-checkout for MOM.res.nc")
  endif()
  execute_process(
    COMMAND git -C "${_SOCA_CLONE_DIR}" lfs pull --include=test/Data/72x35x25/restarts/MOM.res.nc
    RESULT_VARIABLE _lfs_rc
    OUTPUT_QUIET)
  if(NOT _lfs_rc EQUAL 0)
    message(FATAL_ERROR "Failed to pull LFS file MOM.res.nc")
  endif()
  file(COPY "${_SOCA_CLONE_DIR}/test/Data/72x35x25/restarts/MOM.res.nc"
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/data_generated/mom6")
  file(REMOVE_RECURSE "${_SOCA_CLONE_DIR}")
  message(STATUS "Fetched MOM.res.nc -> ${_MOM_RES_NC}")
else()
  message(STATUS "MOM.res.nc already exists, skipping fetch")
endif()

# Test data generation fixture
add_test(NAME ijedi_setup_mom6_testdata
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/data/mom6/gen_mom6_testdata.py
        ${CMAKE_CURRENT_BINARY_DIR}/data_generated/mom6)

add_test(NAME ijedi_setup_obs_testdata
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/data/obs/gen_obs_testdata.py
        ${CMAKE_CURRENT_BINARY_DIR}/data_generated/obs)

# Unit tests
foreach(TIER RANGE 1 3)
    # Find all YAML files for this tier
    file(GLOB TIER_YAML_FILES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier${TIER}_unit*.yaml")

    # Print YAML FILES
    list(LENGTH TIER_YAML_FILES NUM_FILES)
    message(STATUS "Found ${NUM_FILES} YAML files for tier ${TIER}:")

    foreach(YAML_FILE ${TIER_YAML_FILES})
        # Get the filename without directory and extension
        get_filename_component(TEST_NAME ${YAML_FILE} NAME_WE)

        # Extract components: tier<N>_<type>_<name>
        string(REPLACE "_" ";" NAME_PARTS ${TEST_NAME})
        list(LENGTH NAME_PARTS NUM_PARTS)

        # Extract executalbe name from third element of NAME_PARTS. Executalb is first element when split by -
        list(GET NAME_PARTS 2 EXECUTABLE_NAME)
        string(REPLACE "-" ";" EXECUTABLE_PARTS ${EXECUTABLE_NAME})
        list(GET EXECUTABLE_PARTS 0 EXECUTABLE_BASE)

        # Ensure we have exactly 3 parts (tier, type, name)
        if(${NUM_PARTS} EQUAL 3)
            # Use absolute YAML path so the test works regardless of WORKING_DIRECTORY
            add_test(NAME ijedi_${TEST_NAME}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6 ${IJEDI_MPI_ARGS}
                $<TARGET_FILE:ijedi_test.x> ${EXECUTABLE_BASE} ${YAML_FILE})

            message(STATUS "  Added tier ${TIER} test: ijedi_${TEST_NAME}")
        else()
            # Fail with useful message
            message(FATAL_ERROR "Test file ${YAML_FILE} does not follow naming convention tier<N>_<type>_<name>.yaml")
        endif()
    endforeach()
endforeach()

# MOM6 structured-grid check against SOCA reference gridspec
add_test(NAME ijedi_check_mom6_structured_grid
    COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/data/mom6/check_mom6_structured_grid.py
        ${CMAKE_CURRENT_BINARY_DIR}/mom6_structured_grid.nc
        ${CMAKE_CURRENT_BINARY_DIR}/data_generated/mom6/soca_gridspec.72x35x25.nc)

# HofX3D tests
add_test(NAME ijedi_hofx3d_mom6
         COMMAND $<TARGET_FILE:ijedi.x> hofx3d
                 ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_hofx3d-mom6.yaml)

add_test(NAME ijedi_hofx3d_fv3
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6 ${IJEDI_MPI_ARGS}
                 $<TARGET_FILE:ijedi.x> hofx3d
                 ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_hofx3d-fv3.yaml)

# Comment out until https://github.com/JCSDA-internal/oops/pull/3255 is merged
#add_test(NAME ijedi_hofx3d_coupled
#         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6 ${IJEDI_MPI_ARGS}
#                 $<TARGET_FILE:ijedi_coupled.x> hofx3d
#                 ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_hofx3d-coupled.yaml)

# Variational tests
add_test(NAME ijedi_3dvar_mom6
         COMMAND $<TARGET_FILE:ijedi.x> var
                 ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_3dvar-mom6.yaml)

# Convert state tests
add_test(NAME ijedi_convertstate-mom6_to_mom6
	 COMMAND $<TARGET_FILE:ijedi.x> convertstate
                 ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_convertstate-mom6_to_mom6.yaml)
add_test(NAME ijedi_convertstate-mom6_to_fv3
	 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6 ${IJEDI_MPI_ARGS}
	         $<TARGET_FILE:ijedi.x> convertstate
                 ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_convertstate-mom6_to_fv3.yaml)


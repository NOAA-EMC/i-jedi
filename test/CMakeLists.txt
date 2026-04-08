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

        # Skip this here since we have a special test that downloads the mpas test data and we do not want to double register 
        if(TEST_NAME STREQUAL "tier1_unit_geometry-mpas")
            continue()
        endif()

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


# MPAS geometry test
# ----------------------------------

set(MPAS_TEST_WORKDIR ${CMAKE_CURRENT_BINARY_DIR}/testrun/geometry-mpas)
set(MPAS_DATA_CLONE_DIR ${CMAKE_CURRENT_BINARY_DIR}/testrun/geometry-mpas)
set(MPAS_GEOMETRY_YAML ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-mpas.yaml)
file(MAKE_DIRECTORY ${MPAS_TEST_WORKDIR})
file(CREATE_LINK
    "${CMAKE_CURRENT_BINARY_DIR}/data"
    "${MPAS_TEST_WORKDIR}/data"
    SYMBOLIC
)
# Fetch MPAS test data as a setup fixture test so ctest controls execution order.
add_test(NAME ijedi_setup_mpas_testdata
    COMMAND ${CMAKE_COMMAND}
        -DCLONE_DIR=${MPAS_DATA_CLONE_DIR}
        -DDEST_DIR=${MPAS_TEST_WORKDIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/fetch_mpas_testdata.cmake)
set_tests_properties(ijedi_setup_mpas_testdata PROPERTIES
    FIXTURES_SETUP mpas_testdata)

# Run geometry test
add_test(NAME ijedi_tier1_unit_geometry-mpas
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6 ${IJEDI_MPI_ARGS}
    $<TARGET_FILE:ijedi_test.x> ${EXECUTABLE_BASE} ${YAML_FILE})
set_tests_properties(ijedi_tier1_unit_geometry-mpas PROPERTIES
    FIXTURES_REQUIRED  mpas_testdata
    WORKING_DIRECTORY  ${MPAS_TEST_WORKDIR}
    FIXTURES_SETUP     mpas_grid_output
    ENVIRONMENT        "OMP_NUM_THREADS=1")

# Integration tests
foreach(TIER RANGE 1 3)
    # Find all YAML files for this tier
    file(GLOB TIER_YAML_FILES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier${TIER}_integration*.yaml")

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
                $<TARGET_FILE:ijedi.x> ${EXECUTABLE_BASE} ${YAML_FILE})

            message(STATUS "  Added tier ${TIER} test: ijedi_${TEST_NAME}")
        else()
            # Fail with useful message
            message(FATAL_ERROR "Test file ${YAML_FILE} does not follow naming convention tier<N>_<type>_<name>.yaml")
        endif()
    endforeach()
endforeach()

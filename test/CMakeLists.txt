# Build the test executable
# -------------------------
add_subdirectory(mains)

# Python3 is needed for test data generation
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Generate MOM6 test data
# -----------------------
# ocean_hgrid.nc and ocean_topog.nc are rebuilt from CDL files at test time
# using ncgen, avoiding git-lfs for binary NetCDF assets.
set(MOM6_TESTDATA_DIR "${CMAKE_CURRENT_BINARY_DIR}/testdata/mom6")
set(MOM6_WORKDIR      "${CMAKE_CURRENT_BINARY_DIR}/workdir/mom6")
file(MAKE_DIRECTORY   "${MOM6_WORKDIR}")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-mom6.yaml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/testinput/tier1_unit_geometry-mom6.yaml"
    @ONLY)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_state-mom6.yaml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/testinput/tier1_unit_state-mom6.yaml"
    @ONLY)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_increment-mom6.yaml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/testinput/tier1_unit_increment-mom6.yaml"
    @ONLY)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_errorcovariance-mom6.yaml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/testinput/tier1_unit_errorcovariance-mom6.yaml"
    @ONLY)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_dirac-mom6.yaml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/testinput/tier1_unit_dirac-mom6.yaml"
    @ONLY)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_dirac-diffusion-params-mom6.yaml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/testinput/tier1_unit_dirac-diffusion-params-mom6.yaml"
    @ONLY)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_dirac-diffusion-mom6.yaml.in"
    "${CMAKE_CURRENT_BINARY_DIR}/testinput/tier1_unit_dirac-diffusion-mom6.yaml"
    @ONLY)
add_custom_command(
    OUTPUT "${MOM6_TESTDATA_DIR}/INPUT/ocean_hgrid.nc"
           "${MOM6_TESTDATA_DIR}/INPUT/ocean_topog.nc"
           "${MOM6_TESTDATA_DIR}/soca_gridspec.72x35x25.nc"
    COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/testdata/mom6/gen_mom6_testdata.py"
            "${MOM6_TESTDATA_DIR}"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/testdata/mom6/gen_mom6_testdata.py"
            "${CMAKE_CURRENT_SOURCE_DIR}/testdata/mom6/ocean_hgrid.cdl"
            "${CMAKE_CURRENT_SOURCE_DIR}/testdata/mom6/ocean_topog.cdl"
            "${CMAKE_CURRENT_SOURCE_DIR}/testdata/mom6/soca_gridspec.72x35x25.cdl"
    COMMENT "Generating MOM6 geometry test data from CDL"
    VERBATIM)
add_custom_target(mom6_testdata ALL
    DEPENDS "${MOM6_TESTDATA_DIR}/INPUT/ocean_hgrid.nc")

# CTest fixture so MOM6 tests wait for data generation
add_test(NAME mom6_testdata_setup
    COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/testdata/mom6/gen_mom6_testdata.py"
            "${MOM6_TESTDATA_DIR}")
set_tests_properties(mom6_testdata_setup PROPERTIES
    FIXTURES_SETUP mom6_testdata_fixture)

# Automatically discover and add tests organized by tier
# -------------------------------------------------------
# All tests are built; use ctest -R tier<N> to run specific tiers
# Note: CONFIGURE_DEPENDS tells CMake to recheck for new/removed test files
# Extra MPI arguments, e.g. "--map-by :OVERSUBSCRIBE" for oversubscription.
# Cached so it persists across reconfigures; initialised from env var on first configure.
set(IJEDI_MPI_ARGS "$ENV{IJEDI_MPI_ARGS}" CACHE STRING
    "Extra mpiexec arguments (e.g. --map-by :OVERSUBSCRIBE)")

if(IJEDI_MPI_ARGS)
    separate_arguments(IJEDI_MPI_ARGS)
endif()

foreach(TIER RANGE 1 3)
    # Find all YAML files for this tier (source tree + configured files in build tree)
    file(GLOB TIER_YAML_FILES CONFIGURE_DEPENDS
         "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier${TIER}_*.yaml"
         "${CMAKE_CURRENT_BINARY_DIR}/testinput/tier${TIER}_*.yaml")

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
            # Get relative path for ARGS
            file(RELATIVE_PATH YAML_REL_PATH ${CMAKE_CURRENT_BINARY_DIR} ${YAML_FILE})

            # Determine MPI proc count from filename suffix (-NprocP), default 6
            set(NPROC 6)
            if(EXECUTABLE_NAME MATCHES "-([0-9]+)$")
                set(NPROC ${CMAKE_MATCH_1})
            endif()

            # MOM6 tests run in a separate workdir, so use absolute YAML path
            if(TEST_NAME MATCHES "mom6")
                set(YAML_PATH ${YAML_FILE})
            else()
                set(YAML_PATH ${YAML_REL_PATH})
            endif()

            # Add the test
            add_test(NAME ijedi_${TEST_NAME}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${IJEDI_MPI_ARGS}
                $<TARGET_FILE:ijedi_test.x> ${EXECUTABLE_BASE} ${YAML_PATH})

            # MOM6 test requires generated test data and sets up output fixture
            if(TEST_NAME MATCHES "mom6")
                set_tests_properties(ijedi_${TEST_NAME} PROPERTIES
                    FIXTURES_REQUIRED  mom6_testdata_fixture
                    FIXTURES_SETUP     mom6_geometry_fixture
                    WORKING_DIRECTORY  "${MOM6_WORKDIR}")
            endif()

            message(STATUS "  Added tier ${TIER} test: ijedi_${TEST_NAME} (${NPROC} procs)")
        else()
            # Fail with useful message
            message(FATAL_ERROR "Test file ${YAML_FILE} does not follow naming convention tier<N>_<type>_<name>.yaml")
        endif()
    endforeach()
endforeach()

# Diffusion tests: params generation must run before the dirac test.
# These overrides replace the generic MOM6 fixture set from the loop above.
set_tests_properties(ijedi_tier1_unit_dirac-diffusion-params-mom6 PROPERTIES
    FIXTURES_REQUIRED  mom6_testdata_fixture
    FIXTURES_SETUP     "mom6_geometry_fixture;mom6_diffusion_params_fixture"
    WORKING_DIRECTORY  "${MOM6_WORKDIR}")
set_tests_properties(ijedi_tier1_unit_dirac-diffusion-mom6 PROPERTIES
    FIXTURES_REQUIRED  "mom6_testdata_fixture;mom6_diffusion_params_fixture"
    FIXTURES_SETUP     mom6_geometry_fixture
    WORKING_DIRECTORY  "${MOM6_WORKDIR}")

# Structured-grid regression test: compare geometry output against SOCA reference.
# Registered after the foreach loop so its test number follows the geometry test.
add_test(NAME ijedi_tier1_unit_geometry-mom6_check
    COMMAND ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_SOURCE_DIR}/testdata/mom6/check_mom6_structured_grid.py"
            "${MOM6_WORKDIR}/mom6_structured_grid.nc"
            "${MOM6_TESTDATA_DIR}/soca_gridspec.72x35x25.nc")
set_tests_properties(ijedi_tier1_unit_geometry-mom6_check PROPERTIES
    FIXTURES_REQUIRED mom6_geometry_fixture)

# ---------------------------------------------------------------------------
# Tier-2 high-resolution MOM6 geometry test (optional)
# ---------------------------------------------------------------------------
# Requires user-supplied ocean_hgrid.nc + ocean_topog.nc (flat directory,
# NOT in an INPUT/ subdir).  CMake creates INPUT/ symlinks in a staging area
# so the geometry code's expected layout is satisfied without copying data.
#
# Enable at configure time:
#   cmake -DMOM6_HIRES_DIR=/path/to/dir/containing/ocean_hgrid.nc ..
#
# Grid parameters (0.25° OM4, 6 MPI tasks) are hardcoded in the YAML files.
# ---------------------------------------------------------------------------
set(MOM6_HIRES_DIR "" CACHE PATH
    "Directory containing ocean_hgrid.nc + ocean_topog.nc for tier-2 MOM6 test")

if(MOM6_HIRES_DIR)
    # Hardcoded: 3x2 layout = 6 MPI tasks
    set(MOM6_HIRES_NPROC 6)

    # Staging area: geometry expects INPUT/ subdir; create it with symlinks
    set(MOM6_HIRES_STAGING "${CMAKE_CURRENT_BINARY_DIR}/testdata/mom6_hires")
    file(MAKE_DIRECTORY "${MOM6_HIRES_STAGING}/INPUT")
    foreach(_grid_file ocean_hgrid.nc ocean_topog.nc)
        if(NOT EXISTS "${MOM6_HIRES_STAGING}/INPUT/${_grid_file}")
            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${MOM6_HIRES_DIR}/${_grid_file}"
                "${MOM6_HIRES_STAGING}/INPUT/${_grid_file}")
        endif()
    endforeach()

    # Configure YAML into hires/ subdir (not touched by the auto-discovery loop)
    # The only substituted variable is @MOM6_HIRES_STAGING@ (the data path).
    set(MOM6_HIRES_YAML_DIR "${CMAKE_CURRENT_BINARY_DIR}/testinput/hires")
    set(MOM6_HIRES_YAML     "${MOM6_HIRES_YAML_DIR}/tier2_unit_geometry-mom6.yaml")
    set(MOM6_HIRES_WORKDIR  "${CMAKE_CURRENT_BINARY_DIR}/workdir/mom6_hires")
    file(MAKE_DIRECTORY "${MOM6_HIRES_YAML_DIR}" "${MOM6_HIRES_WORKDIR}")

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier2_unit_geometry-mom6.yaml.in"
        "${MOM6_HIRES_YAML}"
        @ONLY)

    add_test(NAME ijedi_tier2_unit_geometry-mom6
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MOM6_HIRES_NPROC}
                ${IJEDI_MPI_ARGS}
                $<TARGET_FILE:ijedi_test.x> geometry "${MOM6_HIRES_YAML}")
    set_tests_properties(ijedi_tier2_unit_geometry-mom6 PROPERTIES
        WORKING_DIRECTORY  "${MOM6_HIRES_WORKDIR}"
        FIXTURES_SETUP     mom6_hires_geometry_fixture)

    # Optional SOCA reference gridspec for the structured-grid regression check
    set(MOM6_HIRES_SOCA_GRIDSPEC "${MOM6_HIRES_DIR}/soca_gridspec.nc" CACHE PATH
        "SOCA reference gridspec for tier-2 structured-grid check (leave empty to skip)")

    if(MOM6_HIRES_SOCA_GRIDSPEC AND EXISTS "${MOM6_HIRES_SOCA_GRIDSPEC}")
        # areaT tolerance is relaxed: SOCA's area comes from MOM6's runtime spherical-geometry
        # initialization, while i-jedi sums supergrid sub-cells → systematic ~0.8% difference.
        add_test(NAME ijedi_tier2_unit_geometry-mom6_check
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/testdata/mom6/check_mom6_structured_grid.py"
                    "${MOM6_HIRES_WORKDIR}/mom6_hires_structured_grid.nc"
                    "${MOM6_HIRES_SOCA_GRIDSPEC}"
                    "--rtol-area=1e-2")
        set_tests_properties(ijedi_tier2_unit_geometry-mom6_check PROPERTIES
            FIXTURES_REQUIRED mom6_hires_geometry_fixture)
        message(STATUS "  soca ref  : ${MOM6_HIRES_SOCA_GRIDSPEC}")
    else()
        message(STATUS "  soca ref  : not found — structured-grid check skipped")
        message(STATUS "  set -DMOM6_HIRES_SOCA_GRIDSPEC=/path/to/soca_gridspec.nc to enable")
    endif()

    # ---------------------------------------------------------------------------
    # Tier-2 high-resolution MOM6 diffusion tests
    # ---------------------------------------------------------------------------
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier2_unit_dirac-diffusion-params-mom6.yaml.in"
        "${MOM6_HIRES_YAML_DIR}/tier2_unit_dirac-diffusion-params-mom6.yaml"
        @ONLY)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier2_unit_dirac-diffusion-mom6.yaml.in"
        "${MOM6_HIRES_YAML_DIR}/tier2_unit_dirac-diffusion-mom6.yaml"
        @ONLY)

    add_test(NAME ijedi_tier2_unit_dirac-diffusion-params-mom6
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MOM6_HIRES_NPROC}
                ${IJEDI_MPI_ARGS}
                $<TARGET_FILE:ijedi_test.x> dirac
                "${MOM6_HIRES_YAML_DIR}/tier2_unit_dirac-diffusion-params-mom6.yaml")
    set_tests_properties(ijedi_tier2_unit_dirac-diffusion-params-mom6 PROPERTIES
        WORKING_DIRECTORY "${MOM6_HIRES_WORKDIR}"
        FIXTURES_SETUP    mom6_hires_diffusion_params_fixture)

    add_test(NAME ijedi_tier2_unit_dirac-diffusion-mom6
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MOM6_HIRES_NPROC}
                ${IJEDI_MPI_ARGS}
                $<TARGET_FILE:ijedi_test.x> dirac
                "${MOM6_HIRES_YAML_DIR}/tier2_unit_dirac-diffusion-mom6.yaml")
    set_tests_properties(ijedi_tier2_unit_dirac-diffusion-mom6 PROPERTIES
        WORKING_DIRECTORY "${MOM6_HIRES_WORKDIR}"
        FIXTURES_REQUIRED mom6_hires_diffusion_params_fixture)

    message(STATUS "Tier-2 high-res MOM6 tests: ENABLED (6 MPI tasks, 0.25° OM4)")
    message(STATUS "  data: ${MOM6_HIRES_DIR}")
    message(STATUS "  run -R tier2 to execute")
else()
    message(STATUS "Tier-2 high-res MOM6 tests: DISABLED")
    message(STATUS "  set -DMOM6_HIRES_DIR=/path/to/dir/with/ocean_hgrid.nc to enable")
endif()

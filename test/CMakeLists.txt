# Build the test executable
# -------------------------
add_subdirectory(mains)

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

# Test data generation
ecbuild_add_test(TARGET  ijedi_setup_mom6_testdata
                 TYPE    SCRIPT
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/gen_mom6_testdata.py
                 ARGS    ${CMAKE_CURRENT_SOURCE_DIR}/data/mom6 ${CMAKE_CURRENT_BINARY_DIR}/data_generated/mom6)

ecbuild_add_test(TARGET  ijedi_setup_obs_testdata
                 TYPE    SCRIPT
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/gen_obs_testdata.py
                 ARGS    ${CMAKE_CURRENT_SOURCE_DIR}/data/obs ${CMAKE_CURRENT_BINARY_DIR}/data_generated/obs)

ecbuild_add_test(TARGET ijedi_setup_fv3_testdata
                 TYPE SCRIPT
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/gen_fv3_testdata.py
                 ARGS    ${CMAKE_CURRENT_BINARY_DIR}/data_generated/fv3)

# Geometry tests
ecbuild_add_test( TARGET  ijedi_geometry_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometry ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-fv3.yaml)

ecbuild_add_test( TARGET  ijedi_geometry_mom6
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometry ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-mom6.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata)

# MOM6 structured-grid check against SOCA reference gridspec
ecbuild_add_test(TARGET  ijedi_check_mom6_structured_grid
                 TYPE    SCRIPT
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_mom6_structured_grid.py
                 ARGS    ${CMAKE_CURRENT_BINARY_DIR}/mom6_structured_grid.nc
                         ${CMAKE_CURRENT_BINARY_DIR}/data_generated/mom6/soca_gridspec.72x35x25.nc
                 TEST_DEPENDS ijedi_geometry_mom6)

# Increment tests
ecbuild_add_test( TARGET  ijedi_increment_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    increment ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_increment-fv3.yaml)

# State tests
ecbuild_add_test( TARGET  ijedi_state_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    state ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_state-fv3.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata)

ecbuild_add_test( TARGET  ijedi_state_mom6
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    state ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_state-mom6.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata)

# H(x) tests
ecbuild_add_test( TARGET  ijedi_hofx3d_mom6
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    hofx3d ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_hofx3d-mom6.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata ijedi_setup_obs_testdata)

ecbuild_add_test( TARGET  ijedi_hofx3d_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    hofx3d ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_hofx3d-fv3.yaml
                  TEST_DEPENDS ijedi_setup_obs_testdata)

ecbuild_add_test( TARGET  ijedi_hofx3d_coupled
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_coupled.x>
                  ARGS    hofx3d ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_hofx3d-coupled.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata ijedi_setup_obs_testdata)

# Variational tests (depends on diffusion operator being calibrated first)
ecbuild_add_test( TARGET  ijedi_diffusion_mom6
                  MPI     2
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    errortoolbox ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_diffusion-mom6.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata)

ecbuild_add_test( TARGET  ijedi_3dvar_mom6
                  MPI     2
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    var ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_3dvar-mom6.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata ijedi_setup_obs_testdata ijedi_diffusion_mom6)

# Convert state tests
ecbuild_add_test( TARGET  ijedi_convertstate_mom6_to_mom6
                  MPI     1
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    convertstate ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_convertstate-mom6_to_mom6.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata)

ecbuild_add_test( TARGET  ijedi_convertstate_mom6_to_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    convertstate ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_convertstate-mom6_to_fv3.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata ijedi_setup_fv3_testdata)

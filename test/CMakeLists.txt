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
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data_generated/fv3)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data_generated/aimodel)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data_generated/neuralgcm)

# ###############################################################################
# NeuralGCM tests
#
# NeuralGCM runs in its own python virtual environment, deliberately outside the JEDI build:
# it needs numpy 2.x where spack-stack provides 1.26, so the two cannot be merged. Create the
# environment with tools/setup_python_model_env.sh and point IJEDI_PYTHON_ENV at it, either
# here or in the environment at test time. Without it these tests report as skipped rather
# than failed - a build on a machine with no model environment is not broken.
# ###############################################################################
set(IJEDI_PYTHON_ENV "$ENV{IJEDI_PYTHON_ENV}" CACHE PATH
    "Virtual environment holding the python forecast models (NeuralGCM and friends)")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tools/run_if_neuralgcm.sh
               ${CMAKE_CURRENT_BINARY_DIR}/run_if_neuralgcm.sh @ONLY)

if(NOT IJEDI_PYTHON_ENV)
  message(STATUS "IJEDI_PYTHON_ENV is not set: the NeuralGCM tests will report as skipped")
else()
  message(STATUS "NeuralGCM python environment: ${IJEDI_PYTHON_ENV}")
endif()

# Link the directory 'testref' to the test directory of the build
# ---------------------------------------------------------------
file(CREATE_LINK
    "${CMAKE_CURRENT_SOURCE_DIR}/testref"
    "${CMAKE_CURRENT_BINARY_DIR}/testref"
    SYMBOLIC
)

# Test data generation
ecbuild_add_test(TARGET  ijedi_setup_mom6_testdata
                 TYPE    SCRIPT
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/gen_mom6_testdata.py
                 ARGS    ${CMAKE_CURRENT_SOURCE_DIR}/data/mom6 ${CMAKE_CURRENT_BINARY_DIR}/data_generated/mom6)

ecbuild_add_test(TARGET  ijedi_setup_obs_testdata
                 TYPE    SCRIPT
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/gen_obs_testdata.py
                 ARGS    ${CMAKE_CURRENT_SOURCE_DIR}/data/obs ${CMAKE_CURRENT_BINARY_DIR}/data_generated/obs)

ecbuild_add_test(TARGET  ijedi_setup_aimodel_testdata
                 TYPE    SCRIPT
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/gen_aimodel_testdata.py
                 ARGS    ${CMAKE_CURRENT_BINARY_DIR}/data_generated/aimodel)

ecbuild_add_test(TARGET ijedi_setup_fv3_testdata
                 TYPE SCRIPT
                 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tools/gen_fv3_testdata.py
                 ARGS    ${CMAKE_CURRENT_BINARY_DIR}/data_generated/fv3)

ecbuild_add_test(TARGET  ijedi_setup_neuralgcm_testdata
                 TYPE    SCRIPT
                 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/run_if_neuralgcm.sh
                 ARGS    PYTHON
                         ${CMAKE_CURRENT_SOURCE_DIR}/../tools/stage_neuralgcm_data.py
                         ${CMAKE_CURRENT_BINARY_DIR}/data_generated/neuralgcm
                         --hours 6 --seed 0)
set_tests_properties(ijedi_setup_neuralgcm_testdata PROPERTIES SKIP_RETURN_CODE 77)

# Geometry tests
ecbuild_add_test( TARGET  ijedi_geometry_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometry ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-fv3.yaml)

ecbuild_add_test( TARGET  ijedi_geometry_atlas
                  MPI     2
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometry ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-atlas.yaml)

ecbuild_add_test( TARGET  ijedi_geometry_neuralgcm
                  MPI     2
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometry
                          ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-neuralgcm.yaml)

ecbuild_add_test( TARGET  ijedi_geometry_gsibec
                  MPI     2
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometry ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-gsibec.yaml)

ecbuild_add_test( TARGET  ijedi_geometry_mom6
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometry ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-mom6.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata)

# OOPS GeometryIterator interface tests
ecbuild_add_test( TARGET  ijedi_geometryiterator_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometryiterator ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-fv3.yaml)

ecbuild_add_test( TARGET  ijedi_geometryiterator_mom6
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometryiterator ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-mom6.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata)

ecbuild_add_test( TARGET  ijedi_geometryiterator_neuralgcm_3diter
                  MPI     2
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometryiterator
                          ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-neuralgcm.yaml)

ecbuild_add_test( TARGET  ijedi_geometryiterator_fv3_3diter
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometryiterator
                          ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-fv3-3diter.yaml)

ecbuild_add_test( TARGET  ijedi_geometryiterator_mom6_3diter
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    geometryiterator
                          ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_geometry-mom6-3diter.yaml
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
                  TEST_DEPENDS ijedi_setup_mom6_testdata ijedi_geometry_cache_mom6)

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

# Geometry cache test
ecbuild_add_test( TARGET  ijedi_geometry_cache_mom6
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    geometry_cache ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_geometry_cache-mom6.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata)

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
                  TEST_DEPENDS ijedi_setup_mom6_testdata ijedi_setup_obs_testdata ijedi_diffusion_mom6 ijedi_geometry_cache_mom6)

# 3DVar comparison: FV3 increment with interpolation to gaussian (Variant A) vs
# native atlas/gaussian increment (Variant B).
ecbuild_add_test( TARGET  ijedi_3dvar_fv3_interp_gauss
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    var ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_3dvar-fv3_interp_gauss.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata ijedi_setup_obs_testdata)

ecbuild_add_test( TARGET  ijedi_3dvar_fv3_inc_gauss
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    var ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_3dvar-fv3_inc_gauss.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata ijedi_setup_obs_testdata)

# LETKF tests
ecbuild_add_test( TARGET  ijedi_letkf_observer_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    letkf ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_letkf_observer-fv3.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata ijedi_setup_obs_testdata)

ecbuild_add_test( TARGET  ijedi_letkf_solver_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    letkf ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_letkf_solver-fv3.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata ijedi_setup_obs_testdata ijedi_letkf_observer_fv3)

ecbuild_add_test( TARGET  ijedi_letkf_observer_fv3-restart
                  MPI     1
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    letkf ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_letkf_observer-fv3-restart.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata ijedi_setup_obs_testdata)

ecbuild_add_test( TARGET  ijedi_letkf_solver_fv3-restart
                  MPI     1
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    letkf ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_letkf_solver-fv3-restart.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata ijedi_setup_obs_testdata ijedi_letkf_observer_fv3-restart)

# Dirac tests
ecbuild_add_test( TARGET  ijedi_dirac_fv3_interp_gauss
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    errortoolbox ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_dirac-fv3_interp_gauss.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata)

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

# A NeuralGCM forecast, on NeuralGCM's own native grid, checked field by field against a
# rollout that ran entirely in Python and never went through the bridge.
ecbuild_add_test( TARGET  ijedi_forecast_neuralgcm
                  MPI     2
                  OMP     1
                  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/run_if_neuralgcm.sh
                  ARGS    $<TARGET_FILE:ijedi.x> forecast
                          ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_forecast-neuralgcm.yaml
                  TEST_DEPENDS ijedi_setup_neuralgcm_testdata)
set_tests_properties(ijedi_forecast_neuralgcm PROPERTIES SKIP_RETURN_CODE 77)

ecbuild_add_test( TARGET  ijedi_check_neuralgcm_forecast
                  TYPE    SCRIPT
                  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/run_if_neuralgcm.sh
                  ARGS    PYTHON
                          ${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_neuralgcm_forecast.py
                          ${CMAKE_CURRENT_BINARY_DIR}/forecast-neuralgcm-tl63.nc
                          ${CMAKE_CURRENT_BINARY_DIR}/data_generated/neuralgcm/neuralgcm_ref_tl63.nc
                  TEST_DEPENDS ijedi_forecast_neuralgcm)
set_tests_properties(ijedi_check_neuralgcm_forecast PROPERTIES SKIP_RETURN_CODE 77)

# Linear model: the oops battery against i-jedi's persistence linear model. Exercises the
# LinearModel factory, dispatcher, increment handling and the forward and backward clock,
# with no python involved. Persistence is exactly linear and its own adjoint, so both the
# TL and AD cases should hold to machine precision.
ecbuild_add_test( TARGET  ijedi_linearmodel_persistence
                  MPI     2
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi_test.x>
                  ARGS    linearmodel
                          ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_unit_linearmodel-persistence.yaml
                  TEST_DEPENDS ijedi_setup_aimodel_testdata)

# NeuralGCM's tangent linear and adjoint, checked against each other by the dot-product
# identity. Derived by automatic differentiation, so this should hold to rounding.
ecbuild_add_test( TARGET  ijedi_tlad_neuralgcm
                  TYPE    SCRIPT
                  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/run_if_neuralgcm.sh
                  ARGS    PYTHON
                          ${CMAKE_CURRENT_SOURCE_DIR}/tools/test_neuralgcm_tlad.py
                          ${CMAKE_CURRENT_BINARY_DIR}/data_generated/neuralgcm/neuralgcm_ic_tl63.nc
                  TEST_DEPENDS ijedi_setup_neuralgcm_testdata)
set_tests_properties(ijedi_tlad_neuralgcm PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 1800)

# Python-model interface, exercised with the stub backend (no Python required)
ecbuild_add_test( TARGET  ijedi_forecast_python_stub
                  MPI     4
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    forecast
                          ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_forecast-python-stub.yaml
                  TEST_DEPENDS ijedi_setup_aimodel_testdata)

ecbuild_add_test( TARGET  ijedi_convertstate_aimodel
                  MPI     4
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    convertstate
                          ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_convertstate-aimodel.yaml
                  TEST_DEPENDS ijedi_setup_aimodel_testdata)

ecbuild_add_test( TARGET  ijedi_convertstate_fv3_to_atlas
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    convertstate ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_convertstate-fv3_to_atlas.yaml
                  TEST_DEPENDS ijedi_setup_mom6_testdata ijedi_setup_fv3_testdata)

ecbuild_add_test( TARGET  ijedi_convertstate_fv3_to_fv3
                  MPI     6
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    convertstate ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_convertstate-fv3_to_fv3.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata)

ecbuild_add_test( TARGET  ijedi_convertstate_fv3-restart_to_fv3-restart
                  MPI     1
                  OMP     1
                  COMMAND $<TARGET_FILE:ijedi.x>
                  ARGS    convertstate ${CMAKE_CURRENT_SOURCE_DIR}/testinput/tier1_integration_convertstate-fv3-restart_to_fv3-restart.yaml
                  TEST_DEPENDS ijedi_setup_fv3_testdata)

# fetch_mpas_testdata.cmake
# -------------------------
# Run at CMake configure time (via execute_process in test/CMakeLists.txt).
# Sparse-clones testinput_tier_1/480km/bg from mpas-jedi-data (with LFS) and
# stages the files the MPAS geometry unit test needs inside DEST_DIR, which is
# the WORKING_DIRECTORY of the ctest (geometry-mpas/).
#
# Expected layout after staging (paths relative to DEST_DIR):
#   Data/480km/bg/restart.2018-04-15_00.00.00.nc   ← from mpas-jedi-data LFS
#   x1.2562.graph.info.part.6                       ← generated block-partition
#
# Variables set by the caller (-D):
#   CLONE_DIR  – scratch directory for the cloned repo
#   DEST_DIR   – MPAS test working directory (geometry-mpas/)

cmake_minimum_required(VERSION 3.15)

if(NOT DEFINED CLONE_DIR OR NOT DEFINED DEST_DIR)
  message(FATAL_ERROR "fetch_mpas_testdata.cmake requires CLONE_DIR and DEST_DIR")
endif()

set(REPO_URL   "https://github.com/JCSDA-internal/mpas-jedi-data.git")
set(BRANCH     "develop")
set(SPARSE_PATH "testinput_tier_1/480km/bg")

# ------------------------------------------------------------------
# 1. Sparse-clone (skip if already present)
# ------------------------------------------------------------------
if(NOT EXISTS "${CLONE_DIR}/.git")
  message(STATUS "Cloning mpas-jedi-data (sparse, no blobs)...")
  execute_process(
    COMMAND git clone
      --depth 1 --branch ${BRANCH}
      --filter=blob:none --no-checkout
      ${REPO_URL} ${CLONE_DIR}
    RESULT_VARIABLE rc
  )
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "git clone of mpas-jedi-data failed (rc=${rc})")
  endif()

  execute_process(
    COMMAND git -C ${CLONE_DIR} sparse-checkout init --cone
    RESULT_VARIABLE rc
  )
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "git sparse-checkout init failed (rc=${rc})")
  endif()

  execute_process(
    COMMAND git -C ${CLONE_DIR} sparse-checkout set ${SPARSE_PATH}
    RESULT_VARIABLE rc
  )
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "git sparse-checkout set failed (rc=${rc})")
  endif()

  execute_process(
    COMMAND git -C ${CLONE_DIR} checkout ${BRANCH}
    RESULT_VARIABLE rc
  )
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "git checkout failed (rc=${rc})")
  endif()

  message(STATUS "Fetching LFS objects for ${SPARSE_PATH}...")
  execute_process(
    COMMAND git -C ${CLONE_DIR} lfs pull --include="${SPARSE_PATH}/**"
    RESULT_VARIABLE rc
  )
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "git lfs pull failed (rc=${rc})")
  endif()
else()
  message(STATUS "mpas-jedi-data clone already present at ${CLONE_DIR}, skipping clone.")
endif()

# ------------------------------------------------------------------
# 2. Stage files into the test working directory.
#    Use separate parallel lists and a counter-based loop to avoid
#    CMake's semicolon-as-separator ambiguity inside foreach pairs.
# ------------------------------------------------------------------
set(_SRCS
  "${CLONE_DIR}/testinput_tier_1/480km/bg/restart.2018-04-15_00.00.00.nc"
)
set(_DSTS
  "${DEST_DIR}/Data/480km/bg/restart.2018-04-15_00.00.00.nc"
)

list(LENGTH _SRCS _N)
math(EXPR _LAST "${_N} - 1")
foreach(_i RANGE ${_LAST})
  list(GET _SRCS ${_i} _SRC)
  list(GET _DSTS ${_i} _DST)

  if(NOT EXISTS "${_SRC}")
    message(FATAL_ERROR "Expected file not found after clone: ${_SRC}")
  endif()

  get_filename_component(_DST_DIR "${_DST}" DIRECTORY)
  file(MAKE_DIRECTORY "${_DST_DIR}")

  if(NOT EXISTS "${_DST}")
    file(CREATE_LINK "${_SRC}" "${_DST}" SYMBOLIC)
    message(STATUS "Staged (symlink): ${_DST}")
  else()
    message(STATUS "Already staged: ${_DST}")
  endif()
endforeach()

# ------------------------------------------------------------------
# 3. Generate x1.2562.graph.info.part.6 if absent.
#    MPAS requires a graph partition file when nproc > 1.
#    A simple block partition (427 cells per rank, 6 ranks) is valid.
#    2562 = 6 * 427, so this divides evenly.
# ------------------------------------------------------------------
set(_PART_FILE "${DEST_DIR}/x1.2562.graph.info.part.6")
if(NOT EXISTS "${_PART_FILE}")
  message(STATUS "Generating block partition: ${_PART_FILE}")
  set(_content "")
  set(_NCELLS 2562)
  set(_NPARTS 6)
  set(_PER_PART 427)   # 2562 / 6
  foreach(_cell RANGE 1 ${_NCELLS})
    math(EXPR _part "( ${_cell} - 1 ) / ${_PER_PART}")
    # Clamp to [0, NPARTS-1] for any rounding edge-cases
    if(_part GREATER_EQUAL ${_NPARTS})
      math(EXPR _part "${_NPARTS} - 1")
    endif()
    string(APPEND _content "${_part}\n")
  endforeach()
  file(WRITE "${_PART_FILE}" "${_content}")
else()
  message(STATUS "Already present: ${_PART_FILE}")
endif()

message(STATUS "MPAS test data ready in ${DEST_DIR}")

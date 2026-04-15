#!/bin/bash
# fix_circular_dependency.sh
#
# Breaks the oops <-> vader circular dependency in installed cmake files.
#
# Problem: oops-targets.cmake references "vader" target, and
# vader-targets.cmake references "oops" target. Each validates that
# referenced targets exist at include-time, so neither can load first.
#
# Solution: Patch both *-config.cmake files to directly include the
# OTHER package's targets file right before including their own.
# The targets files have multiple-inclusion guards ("_cmake_targets_defined"),
# so including them twice is safe.
#
# Usage:
#   ./fix_circular_dependency.sh <install_prefix>

set -euo pipefail

INSTALL_DIR="${1:?Usage: $0 <install_prefix>}"

OOPS_CMAKE_DIR="${INSTALL_DIR}/lib/cmake/oops"
VADER_CMAKE_DIR="${INSTALL_DIR}/lib/cmake/vader"
OOPS_CONFIG="${OOPS_CMAKE_DIR}/oops-config.cmake"
VADER_CONFIG="${VADER_CMAKE_DIR}/vader-config.cmake"
OOPS_TARGETS="${OOPS_CMAKE_DIR}/oops-targets.cmake"
VADER_TARGETS="${VADER_CMAKE_DIR}/vader-targets.cmake"

# --- Validate all files exist ---
for f in "$OOPS_CONFIG" "$VADER_CONFIG" "$OOPS_TARGETS" "$VADER_TARGETS"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: File not found: $f" >&2
    exit 1
  fi
done

# --- Remove the cross-package target validation from both targets files ---
# oops-targets.cmake has a block like:
#   foreach(_target "vader" )
#     if(NOT TARGET "${_target}")
#       ...error...
# We comment out those blocks so they don't fail on missing targets.

echo "Patching $OOPS_TARGETS: removing cross-package target check for vader..."
sed -i '/foreach(_target "vader"/,/endforeach/s/^/#PATCHED /' "$OOPS_TARGETS"

echo "Patching $VADER_TARGETS: removing cross-package target check for oops..."
sed -i '/foreach(_target "oops"/,/endforeach/s/^/#PATCHED /' "$VADER_TARGETS"

# --- Patch oops-config.cmake: find vader before loading oops targets ---
# Insert find_package(vader) right before the targets file is loaded.
# We use find_package (not find_dependency) and QUIET so it doesn't fail.
# We also need a re-entrance guard so oops doesn't infinitely recurse.

if ! grep -q '_oops_config_guard' "$OOPS_CONFIG"; then
  echo "Patching $OOPS_CONFIG..."

  # Add guard at the very top
  sed -i '1i\
# Re-entrance guard for circular dependency\
if(_oops_config_guard)\
  return()\
endif()\
set(_oops_config_guard TRUE)' "$OOPS_CONFIG"

  # Add find_package(vader) before the targets file inclusion block
  sed -i '/### insert definitions for IMPORTED targets/i\
# --- Circular dependency fix: ensure vader is loaded before oops targets ---\
if(NOT TARGET vader)\
  find_package(vader QUIET CONFIG)\
endif()\
# --- End circular dependency fix ---\
' "$OOPS_CONFIG"
else
  echo "Skipping $OOPS_CONFIG: already patched."
fi

# --- Patch vader-config.cmake: find oops before loading vader targets ---
if ! grep -q '_vader_config_guard' "$VADER_CONFIG"; then
  echo "Patching $VADER_CONFIG..."

  # Add guard at the very top
  sed -i '1i\
# Re-entrance guard for circular dependency\
if(_vader_config_guard)\
  return()\
endif()\
set(_vader_config_guard TRUE)' "$VADER_CONFIG"

  # Add find_package(oops) before the targets file inclusion block
  sed -i '/### insert definitions for IMPORTED targets/i\
# --- Circular dependency fix: ensure oops is loaded before vader targets ---\
if(NOT TARGET oops)\
  find_package(oops QUIET CONFIG)\
endif()\
# --- End circular dependency fix ---\
' "$VADER_CONFIG"
else
  echo "Skipping $VADER_CONFIG: already patched."
fi

echo ""
echo "Done. Patched files:"
echo "  $OOPS_CONFIG"
echo "  $VADER_CONFIG"
echo "  $OOPS_TARGETS"
echo "  $VADER_TARGETS"
echo ""
echo "Verification:"
grep -n 'PATCHED\|_config_guard\|Circular dependency' \
  "$OOPS_CONFIG" "$VADER_CONFIG" "$OOPS_TARGETS" "$VADER_TARGETS"

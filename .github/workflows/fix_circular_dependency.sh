#!/bin/bash
# fix_circular_dependency.sh
#
# Breaks the oops <-> vader circular dependency in installed cmake files,
# and fixes missing mist headers by patching the oops include paths.
#
# Usage:
#   ./fix_circular_dependency.sh <install_prefix> [oops_source_dir]

set -euo pipefail

INSTALL_DIR="${1:?Usage: $0 <install_prefix> [oops_source_dir]}"
OOPS_SOURCE_DIR="${2:-}"

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

# =====================================================================
# PART 1: Fix oops <-> vader circular dependency
# =====================================================================

echo "=== Part 1: Fixing oops <-> vader circular dependency ==="

# --- Remove the cross-package target validation from both targets files ---
echo "Patching $OOPS_TARGETS: removing cross-package target check for vader..."
sed -i '/foreach(_target "vader"/,/endforeach/s/^/#PATCHED /' "$OOPS_TARGETS"

echo "Patching $VADER_TARGETS: removing cross-package target check for oops..."
sed -i '/foreach(_target "oops"/,/endforeach/s/^/#PATCHED /' "$VADER_TARGETS"

# --- Patch oops-config.cmake ---
if ! grep -q '_oops_config_guard' "$OOPS_CONFIG"; then
  echo "Patching $OOPS_CONFIG..."
  sed -i '1i\
# Re-entrance guard for circular dependency\
if(_oops_config_guard)\
  return()\
endif()\
set(_oops_config_guard TRUE)' "$OOPS_CONFIG"

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

# --- Patch vader-config.cmake ---
if ! grep -q '_vader_config_guard' "$VADER_CONFIG"; then
  echo "Patching $VADER_CONFIG..."
  sed -i '1i\
# Re-entrance guard for circular dependency\
if(_vader_config_guard)\
  return()\
endif()\
set(_vader_config_guard TRUE)' "$VADER_CONFIG"

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

# =====================================================================
# PART 2: Fix missing mist headers
# =====================================================================

echo ""
echo "=== Part 2: Fixing missing mist headers ==="

if [[ -n "$OOPS_SOURCE_DIR" ]]; then
  MIST_DIR="${OOPS_SOURCE_DIR}/mist"
  if [[ -d "$MIST_DIR" ]]; then
    # Copy mist headers into the install include directory
    # mist expects to be found as "mist/base/Geometry.h" etc.
    echo "Installing mist headers from $MIST_DIR to ${INSTALL_DIR}/include/..."
    mkdir -p "${INSTALL_DIR}/include/mist"
    # Copy all header files preserving directory structure
    cd "$OOPS_SOURCE_DIR"
    find mist -name "*.h" -o -name "*.H" -o -name "*.hpp" | while read -r header; do
      dest="${INSTALL_DIR}/include/${header}"
      mkdir -p "$(dirname "$dest")"
      cp "$header" "$dest"
      echo "  Installed: $header"
    done
    cd - > /dev/null
    echo "mist headers installed."
  else
    echo "WARNING: mist directory not found at $MIST_DIR"
    echo "  mist headers were NOT installed. Build may fail."
  fi
else
  echo "No oops source dir provided. Skipping mist header installation."
  echo "  To install mist headers, run:"
  echo "    $0 <install_prefix> <oops_source_dir>"
fi

echo ""
echo "Done."

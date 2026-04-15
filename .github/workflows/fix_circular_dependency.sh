#!/bin/bash
# fix_circular_dependency.sh
#
# Patches oops-import.cmake.in and vader-import.cmake.in to add
# re-entrance guards that break the oops <-> vader circular dependency.
#
# Usage:
#   ./fix_circular_dependency.sh <oops_source_dir> <vader_source_dir>
#
# Example:
#   ./fix_circular_dependency.sh ./oops ./vader

set -euo pipefail

OOPS_DIR="${1:?Usage: $0 <oops_source_dir> <vader_source_dir>}"
VADER_DIR="${2:?Usage: $0 <oops_source_dir> <vader_source_dir>}"

OOPS_IMPORT="${OOPS_DIR}/oops-import.cmake.in"
VADER_IMPORT="${VADER_DIR}/vader-import.cmake.in"

# --- Validate files exist ---
for f in "$OOPS_IMPORT" "$VADER_IMPORT"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: File not found: $f" >&2
    exit 1
  fi
done

# --- Patch oops-import.cmake.in ---

# 1. Add re-entrance guard at the top (after the first comment line)
if ! grep -q '_oops_import_guard' "$OOPS_IMPORT"; then
  echo "Patching $OOPS_IMPORT: adding re-entrance guard..."
  sed -i '1,/^include(CMakeFindDependencyMacro)/{
    /^include(CMakeFindDependencyMacro)/i\
# Re-entrance guard to break circular dependency with vader\
if(_oops_import_guard)\
    return()\
endif()\
set(_oops_import_guard TRUE)\

  }' "$OOPS_IMPORT"
else
  echo "Skipping $OOPS_IMPORT: guard already present."
fi

# 2. Add find_dependency(vader) before the Fortran compiler check
if ! grep -q 'find_dependency(vader)' "$OOPS_IMPORT"; then
  echo "Patching $OOPS_IMPORT: adding find_dependency(vader)..."
  sed -i '/#Export Fortran compiler version/i\
# vader is needed because oops exported targets reference it\
if(NOT vader_FOUND)\
    find_dependency(vader)\
endif()\
' "$OOPS_IMPORT"
else
  echo "Skipping $OOPS_IMPORT: find_dependency(vader) already present."
fi

# --- Patch vader-import.cmake.in ---

# 1. Add re-entrance guard at the top
if ! grep -q '_vader_import_guard' "$VADER_IMPORT"; then
  echo "Patching $VADER_IMPORT: adding re-entrance guard..."
  sed -i '1,/^include(CMakeFindDependencyMacro)/{
    /^include(CMakeFindDependencyMacro)/i\
# Re-entrance guard to break circular dependency with oops\
if(_vader_import_guard)\
    return()\
endif()\
set(_vader_import_guard TRUE)\

  }' "$VADER_IMPORT"
else
  echo "Skipping $VADER_IMPORT: guard already present."
fi

echo "Done. Circular dependency guards applied."

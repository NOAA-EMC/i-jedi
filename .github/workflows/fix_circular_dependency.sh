#!/bin/bash
# fix_circular_dependency.sh
#
# Patches the INSTALLED oops and vader cmake config files to add
# re-entrance guards that break the oops <-> vader circular dependency.
#
# Usage:
#   ./fix_circular_dependency.sh <install_prefix>
#
# Example:
#   ./fix_circular_dependency.sh ./install

set -euo pipefail

INSTALL_DIR="${1:?Usage: $0 <install_prefix>}"

OOPS_IMPORT="${INSTALL_DIR}/lib/cmake/oops/oops-import.cmake"
OOPS_CONFIG="${INSTALL_DIR}/lib/cmake/oops/oops-config.cmake"
VADER_IMPORT="${INSTALL_DIR}/lib/cmake/vader/vader-import.cmake"

# --- Validate files exist ---
for f in "$OOPS_IMPORT" "$VADER_IMPORT" "$OOPS_CONFIG"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: File not found: $f" >&2
    exit 1
  fi
done

# --- Patch vader-import.cmake: add re-entrance guard ---
if ! grep -q '_vader_import_guard' "$VADER_IMPORT"; then
  echo "Patching $VADER_IMPORT: adding re-entrance guard..."
  sed -i '1i\
# Re-entrance guard to break circular dependency with oops\
if(_vader_import_guard)\
    return()\
endif()\
set(_vader_import_guard TRUE)\
' "$VADER_IMPORT"
else
  echo "Skipping $VADER_IMPORT: guard already present."
fi

# --- Patch oops-import.cmake: add re-entrance guard ---
if ! grep -q '_oops_import_guard' "$OOPS_IMPORT"; then
  echo "Patching $OOPS_IMPORT: adding re-entrance guard..."
  sed -i '1i\
# Re-entrance guard to break circular dependency with vader\
if(_oops_import_guard)\
    return()\
endif()\
set(_oops_import_guard TRUE)\
' "$OOPS_IMPORT"
else
  echo "Skipping $OOPS_IMPORT: guard already present."
fi

# --- Patch oops-import.cmake: add find_dependency(vader) ---
if ! grep -q 'find_dependency(vader)' "$OOPS_IMPORT"; then
  echo "Patching $OOPS_IMPORT: adding find_dependency(vader)..."
  # Insert before the Fortran compiler version export block
  sed -i '/#Export Fortran compiler version/i\
# vader is needed because oops exported targets reference it\
if(NOT vader_FOUND)\
    find_dependency(vader)\
endif()\
' "$OOPS_IMPORT"
else
  echo "Skipping $OOPS_IMPORT: find_dependency(vader) already present."
fi

# --- Patch oops-config.cmake: load import file BEFORE targets file ---
# The default ecbuild config loads import, then targets. The problem is
# that oops-targets.cmake references the vader target, and CMake validates
# imported target dependencies when loading the targets file.
# We need vader to be found (via oops-import.cmake) before the targets
# file is included.
#
# Check if the targets file is loaded AFTER the import file (this is the
# default ecbuild layout and should already be correct). If so, our
# find_dependency(vader) in oops-import.cmake will run first.
echo ""
echo "Verifying oops-config.cmake load order..."
IMPORT_LINE=$(grep -n 'IMPORT_FILE\|import' "$OOPS_CONFIG" | head -5)
TARGETS_LINE=$(grep -n 'TARGETS_FILE\|targets' "$OOPS_CONFIG" | head -5)
echo "  Import references:  $IMPORT_LINE"
echo "  Targets references: $TARGETS_LINE"

echo ""
echo "Done. Circular dependency guards applied to installed cmake files."
echo ""
echo "Patched files:"
echo "  $OOPS_IMPORT"
echo "  $VADER_IMPORT"

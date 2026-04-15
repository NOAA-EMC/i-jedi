#!/bin/bash
# fix_circular_dependency.sh
#
# Patches the INSTALLED oops and vader cmake config files to break
# the circular dependency between oops and vader.
#
# The problem: oops-targets.cmake references the vader target, and
# vader-targets.cmake references the oops target. Neither can load
# without the other's targets already defined.
#
# The fix: Before each *-config.cmake loads its own targets file,
# we directly include the OTHER package's targets file (if it exists
# and hasn't been loaded yet). This avoids the full find_package
# re-entry that triggers the cycle.
#
# Usage:
#   ./fix_circular_dependency.sh <install_prefix>

set -euo pipefail

INSTALL_DIR="${1:?Usage: $0 <install_prefix>}"

OOPS_CONFIG="${INSTALL_DIR}/lib/cmake/oops/oops-config.cmake"
VADER_CONFIG="${INSTALL_DIR}/lib/cmake/vader/vader-config.cmake"

# --- Validate files exist ---
for f in "$OOPS_CONFIG" "$VADER_CONFIG"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: File not found: $f" >&2
    exit 1
  fi
done

# --- Find the actual targets files ---
OOPS_TARGETS=$(find "${INSTALL_DIR}/lib/cmake/oops" -name "oops-targets.cmake" -type f | head -1)
VADER_TARGETS=$(find "${INSTALL_DIR}/lib/cmake/vader" -name "vader-targets.cmake" -type f | head -1)

if [[ -z "$OOPS_TARGETS" ]]; then
  echo "ERROR: Could not find oops-targets.cmake under ${INSTALL_DIR}/lib/cmake/oops/" >&2
  exit 1
fi
if [[ -z "$VADER_TARGETS" ]]; then
  echo "ERROR: Could not find vader-targets.cmake under ${INSTALL_DIR}/lib/cmake/vader/" >&2
  exit 1
fi

echo "Found targets files:"
echo "  oops: $OOPS_TARGETS"
echo "  vader: $VADER_TARGETS"

# --- Patch oops-config.cmake ---
# Inject: include vader's targets file directly before oops loads its own
# targets file. This ensures the "vader" imported target exists before
# oops-targets.cmake tries to reference it.
if ! grep -q 'vader-targets.cmake' "$OOPS_CONFIG"; then
  echo "Patching $OOPS_CONFIG..."
  sed -i "/find_file.*oops.*TARGETS_FILE/i\\
# --- Circular dependency fix: pre-load vader targets ---\\
if(NOT TARGET vader)\\
    set(_vader_targets_file \"${VADER_TARGETS}\")\\
    if(EXISTS \"\${_vader_targets_file}\")\\
        include(\"\${_vader_targets_file}\")\\
    endif()\\
    unset(_vader_targets_file)\\
endif()\\
# --- End circular dependency fix ---" "$OOPS_CONFIG"
else
  echo "Skipping $OOPS_CONFIG: vader-targets patch already present."
fi

# --- Patch vader-config.cmake ---
# Same thing in reverse: include oops's targets file before vader loads its own.
if ! grep -q 'oops-targets.cmake' "$VADER_CONFIG"; then
  echo "Patching $VADER_CONFIG..."
  sed -i "/find_file.*vader.*TARGETS_FILE/i\\
# --- Circular dependency fix: pre-load oops targets ---\\
if(NOT TARGET oops)\\
    set(_oops_targets_file \"${OOPS_TARGETS}\")\\
    if(EXISTS \"\${_oops_targets_file}\")\\
        include(\"\${_oops_targets_file}\")\\
    endif()\\
    unset(_oops_targets_file)\\
endif()\\
# --- End circular dependency fix ---" "$VADER_CONFIG"
else
  echo "Skipping $VADER_CONFIG: oops-targets patch already present."
fi

echo ""
echo "Done. Patches applied to:"
echo "  $OOPS_CONFIG"
echo "  $VADER_CONFIG"
echo ""
echo "You can verify with:"
echo "  grep -n 'Circular dependency' $OOPS_CONFIG $VADER_CONFIG"

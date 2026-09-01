#!/usr/bin/env bash
#
# Run a command only if the NeuralGCM Python environment is available, and tell ctest to
# report the test as skipped rather than failed when it is not.
#
# The model lives in its own virtual environment, deliberately outside the JEDI build (see
# tools/setup_python_model_env.sh for why the two cannot be merged). A build on a machine
# where nobody has created that environment is not broken, and these tests should say so
# rather than going red.
#
# Exit code 77 is the conventional "skipped" code and is what SKIP_RETURN_CODE is set to on
# the tests that use this wrapper.
#
# Usage: run_if_neuralgcm.sh <command> [args...]
#   Honours IJEDI_PYTHON_ENV from the environment, falling back to the path baked in at
#   configure time.

set -uo pipefail

SKIP=77

# The default comes from cmake. Note there is deliberately no sentinel comparison against
# the placeholder here: configure_file substitutes every occurrence, so a line testing for
# the un-substituted placeholder would compare the path against itself and always skip.
VENV="${IJEDI_PYTHON_ENV:-@IJEDI_PYTHON_ENV@}"

if [[ -z "${VENV}" ]]; then
  echo "SKIP: no NeuralGCM environment configured."
  echo "      Create one with tools/setup_python_model_env.sh and either re-run cmake with"
  echo "      -DIJEDI_PYTHON_ENV=<prefix> or export IJEDI_PYTHON_ENV=<prefix>."
  exit ${SKIP}
fi

if [[ ! -x "${VENV}/bin/python" ]]; then
  echo "SKIP: ${VENV}/bin/python does not exist."
  exit ${SKIP}
fi

if ! env -u PYTHONPATH "${VENV}/bin/python" -c "import neuralgcm" >/dev/null 2>&1; then
  echo "SKIP: ${VENV} exists but 'import neuralgcm' fails there."
  echo "      Re-run tools/setup_python_model_env.sh ${VENV}"
  exit ${SKIP}
fi

export IJEDI_PYTHON_ENV="${VENV}"

# The venv's interpreter is only known at run time, so a test that needs it writes the
# literal token PYTHON and we substitute it here.
args=()
for arg in "$@"; do
  if [[ "${arg}" == "PYTHON" ]]; then
    args+=("${VENV}/bin/python")
  else
    args+=("${arg}")
  fi
done

exec "${args[@]}"

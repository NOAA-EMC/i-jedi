#!/usr/bin/env bash
#
# Create the Python virtual environment that the i-jedi Python-model adapters run in.
#
# Two things about a spack-stack host make this less obvious than it looks, and both are
# handled here rather than left to the caller:
#
#   1. spack-stack sets PYTHONPATH with dozens of entries. Those land on sys.path ahead of a
#      venv's own site-packages, so "include-system-site-packages = false" isolates nothing.
#      pip then sees spack's packages as already satisfied and silently omits them from the
#      venv, and at run time the venv picks up spack's numpy 1.26 instead of the numpy 2.x
#      that JAX needs. Every python and pip call here clears PYTHONPATH.
#
#   2. pip itself is a spack package reached through PYTHONPATH, so "python3 -m pip" on the
#      base interpreter fails once PYTHONPATH is cleared. Only the venv's own pip, which
#      "python -m venv" bootstraps through ensurepip, is usable.
#
# The venv is built against spack-stack's interpreter rather than the system one. That costs
# nothing and keeps the in-process transport possible later, which needs an exact ABI match.
#
# Usage: setup_python_model_env.sh [prefix]
#   prefix defaults to ${IJEDI_PYTHON_ENV:-$HOME/.local/ijedi/neuralgcm-venv}

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REQUIREMENTS="${HERE}/requirements-neuralgcm.txt"

PREFIX="${1:-${IJEDI_PYTHON_ENV:-$HOME/.local/ijedi/neuralgcm-venv}}"

if [[ ! -f "${REQUIREMENTS}" ]]; then
  echo "ERROR: cannot find ${REQUIREMENTS}" >&2
  exit 1
fi

BASE_PYTHON="$(command -v python3 || true)"
if [[ -z "${BASE_PYTHON}" ]]; then
  echo "ERROR: no python3 on PATH. Load the JEDI environment first." >&2
  exit 1
fi

echo "base interpreter : ${BASE_PYTHON} ($(${BASE_PYTHON} --version 2>&1))"
echo "venv prefix      : ${PREFIX}"

if [[ -x "${PREFIX}/bin/python" ]]; then
  echo "venv already exists; reusing it. Delete ${PREFIX} to rebuild from scratch."
else
  mkdir -p "$(dirname "${PREFIX}")"
  env -u PYTHONPATH "${BASE_PYTHON}" -m venv "${PREFIX}"
fi

env -u PYTHONPATH "${PREFIX}/bin/pip" install --upgrade --quiet pip setuptools wheel
env -u PYTHONPATH "${PREFIX}/bin/pip" install --quiet -r "${REQUIREMENTS}"

echo
echo "verifying..."
env -u PYTHONPATH "${PREFIX}/bin/python" - <<'PYCHECK'
import jax, numpy, neuralgcm
from neuralgcm import demo
model = neuralgcm.PressureLevelModel.from_checkpoint(
    demo.load_checkpoint_tl63_stochastic())
print(f"  neuralgcm {neuralgcm.__version__} | jax {jax.__version__} "
      f"| numpy {numpy.__version__}")
print(f"  checkpoint loads, internal timestep {model.timestep}")
PYCHECK

echo
echo "done. Configure i-jedi against it with:"
echo "  ecbuild -DIJEDI_PYTHON_ENV=${PREFIX} ..."
echo "or export IJEDI_PYTHON_ENV=${PREFIX} before running cmake."

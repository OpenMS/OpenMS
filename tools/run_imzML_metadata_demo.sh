#!/usr/bin/env bash
# Print imzML dataset metadata using pyOpenMS (both reference fixtures by default).
#
# Usage (from anywhere):
#   /path/to/OpenMS/tools/run_imzML_metadata_demo.sh
#   /path/to/OpenMS/tools/run_imzML_metadata_demo.sh /path/to/file.imzML

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_BIN="${ROOT}/OpenMS-build/pyOpenMS/.venv/bin"
PY="${PYOPENMS_PYTHON:-${VENV_BIN}/python3}"
DEMO="${ROOT}/src/pyOpenMS/demo_imzML_metadata.py"

if [[ ! -f "${PY}" ]]; then
  echo "Python not found: ${PY}" >&2
  echo "Set PYOPENMS_PYTHON to your pyOpenMS venv python." >&2
  exit 1
fi

export PYTHONPATH="${ROOT}/OpenMS-build/pyOpenMS${PYTHONPATH:+:${PYTHONPATH}}"
export PYOPENMS_BUILD="${ROOT}/OpenMS-build/pyOpenMS"

# Prefer the pyOpenMS venv over conda/base on PATH (avoids AWS dylib clashes).
exec env PATH="${VENV_BIN}:${PATH}" "${PY}" "${DEMO}" "$@"

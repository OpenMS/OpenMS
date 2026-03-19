#!/usr/bin/env bash

set -euo pipefail

# Resolve repository root from script location.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

is_windows=0
if [[ "${RUNNER_OS:-}" == "Windows" || "${OS:-}" == "Windows_NT" ]]; then
  is_windows=1
fi

if ! command -v pixi >/dev/null 2>&1; then
  if [[ ${is_windows} -eq 1 ]]; then
    if command -v powershell.exe >/dev/null 2>&1; then
      powershell.exe -NoProfile -Command '$ProgressPreference = "SilentlyContinue"; Invoke-RestMethod -Uri https://pixi.sh/install.ps1 | Invoke-Expression'
    else
      echo "Unable to install pixi: powershell.exe not found on Windows runner." >&2
      exit 1
    fi
  else
    curl -fsSL https://pixi.sh/install.sh | bash
  fi

  # pixi installs into ~/.pixi/bin on all runner OSes.
  export PATH="${HOME}/.pixi/bin:${PATH}"
fi

cd "${REPO_ROOT}"
pixi install --frozen

pixi_env="${REPO_ROOT}/.pixi/envs/default"

# Make pixi-managed binaries available in subsequent workflow steps.
if [[ -n "${GITHUB_PATH:-}" ]]; then
  if [[ ${is_windows} -eq 1 ]]; then
    {
      echo "${pixi_env}/bin"
      echo "${pixi_env}/Library/bin"
      echo "${pixi_env}/Library/nsis"
    } >> "${GITHUB_PATH}"
  else
    echo "${pixi_env}/bin" >> "${GITHUB_PATH}"
  fi
fi

if [[ -n "${GITHUB_ENV:-}" ]]; then
  echo "OPENMS_PIXI_BUILD=ON" >> "${GITHUB_ENV}"
  if [[ ${is_windows} -eq 1 ]]; then
    echo "LIB=${pixi_env}/Library/lib${LIB:+;${LIB}}" >> "${GITHUB_ENV}"
  fi
fi

if [[ -n "${GITHUB_OUTPUT:-}" ]]; then
  echo "pixi_env=${pixi_env}" >> "${GITHUB_OUTPUT}"
  echo "openms_pixi_build=ON" >> "${GITHUB_OUTPUT}"

  if [[ ${is_windows} -eq 1 ]]; then
    echo "cmake_prefix=${pixi_env};${pixi_env}/Library;${pixi_env}/Library/lib/cmake;${pixi_env}/Library/lib;${pixi_env}/lib/cmake;${pixi_env}/lib" >> "${GITHUB_OUTPUT}"
  else
    echo "cmake_prefix=${pixi_env};${pixi_env}/lib/cmake;${pixi_env}/share/cmake;${pixi_env}/lib" >> "${GITHUB_OUTPUT}"
  fi
fi

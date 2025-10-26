# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------
## 3.15 to enable new find policies for FindPython (note the discussion in https://gitlab.kitware.com/cmake/cmake/-/issues/21186)
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("pyTOPP" LANGUAGES C CXX)

# Source layout: we copy ONLY the 'pytopp' subpackage to keep 'pyopenms' as a namespace.
set(_pytopp_pkg_src "${CMAKE_CURRENT_SOURCE_DIR}/src/pyopenms/pytopp")

# Where we stage things in the build/install trees
set(_pytopp_site_dir "${CMAKE_CURRENT_BINARY_DIR}/share/pytopp/site")
set(_pytopp_pkg_dst_root "${_pytopp_site_dir}/pyopenms")
file(MAKE_DIRECTORY "${_pytopp_pkg_dst_root}")

# Copy the whole package: share/pytopp/site/pyopenms/pytopp/...
file(COPY "${_pytopp_pkg_src}" DESTINATION "${_pytopp_pkg_dst_root}")

# Write sitecustomize.py so Python auto-imports our bootstrap on startup
file(WRITE "${_pytopp_site_dir}/sitecustomize.py" [=[
# Auto-imported by Python's 'site' module if on sys.path.
# Extends/creates the 'pyopenms' namespace to include the vendored 'pytopp' subtree.

import importlib
import importlib.machinery
import os
import sys
import types

def _ensure_pyopenms_namespace():
    # site root = this file's directory
    site = os.path.dirname(os.path.abspath(__file__))
    pyo_dir = os.path.join(site, "pyopenms")

    # Try to import an installed pyopenms; fall back to creating a stub package.
    try:
        m = importlib.import_module("pyopenms")
    except Exception:
        m = None

    if m is not None:
        # Extend its search path to also include our vendored folder
        try:
            path_list = list(m.__path__)  # may raise if not a pkg
        except Exception:
            path_list = []
        if pyo_dir not in path_list:
            try:
                m.__path__.append(pyo_dir)
            except Exception:
                m.__path__ = path_list + [pyo_dir]
        return

    # No real pyopenms installed: create a minimal namespace package
    m = types.ModuleType("pyopenms")
    m.__package__ = "pyopenms"
    m.__path__ = [pyo_dir]
    m.__spec__ = importlib.machinery.ModuleSpec("pyopenms", loader=None, is_package=True)
    sys.modules["pyopenms"] = m

_ensure_pyopenms_namespace()
]=])

# Path to the shipped CLI inside the vendored site
set(_pytopp_cli_rel "pyopenms/pytopp/tools/pytopp_pyprophet_cli.py")
set(_pytopp_cli_build "${_pytopp_site_dir}/${_pytopp_cli_rel}")

# Make the CLI executable on Unix (so the uv shebang works)
if(UNIX)
  file(CHMOD "${_pytopp_cli_build}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                GROUP_READ GROUP_EXECUTE
                WORLD_READ WORLD_EXECUTE)
endif()

# ---------------------- Wrapper in bin/ (adds PYTHONPATH) -------------------
if(WIN32)
  set(_wrapper "${CMAKE_CURRENT_BINARY_DIR}/bin/pyTOPP-PyProphet.cmd")
  file(WRITE "${_wrapper}" [=[
@echo off
setlocal

REM Resolve wrapper dir
set "WRAPDIR=%~dp0"
if "%WRAPDIR:~-1%"=="\" set "WRAPDIR=%WRAPDIR:~0,-1%"

REM Vendored site root (..\\share\\pytopp\\site)
set "SITE=%WRAPDIR%\..\share\pytopp\site"
set "PYTHONPATH=%SITE%;%PYTHONPATH%"

REM CLI inside vendored package
set "TOOL=%SITE%\pyopenms\pytopp\tools\pytopp_pyprophet_cli.py"

REM Prefer uv explicitly (Windows doesn't honor shebang)
where uv >NUL 2>&1
if %ERRORLEVEL% EQU 0 (
  uv run --script "%TOOL%" %*
  exit /b %ERRORLEVEL%
)

REM Fallback to python with a helpful message
where python >NUL 2>&1
if %ERRORLEVEL% NEQ 0 (
  echo ERROR: Neither "uv" nor "python" found on PATH.
  echo.
  echo Install uv (recommended) and re-run:
  echo   PowerShell:  iwr https://astral.sh/uv/install.ps1 -UseB ^| iex
  echo   winget:      winget install --id AstralSoftware.UV
  echo.
  echo Or install Python and run:
  echo   python "%TOOL%" --help
  exit /b 127
)

python "%TOOL%" %*
]=])
  install(PROGRAMS "${_wrapper}" DESTINATION bin)
else()
  set(_wrapper "${CMAKE_CURRENT_BINARY_DIR}/bin/pyTOPP-PyProphet")
  file(WRITE "${_wrapper}" [=[
#!/usr/bin/env sh
# Resolve wrapper dir
get_script_dir() {
  TARGET_FILE="$0"
  while [ -L "$TARGET_FILE" ]; do
    DIR=$(dirname -- "$TARGET_FILE")
    TARGET_FILE=$(readlink -- "$TARGET_FILE")
    case "$TARGET_FILE" in
      /*) ;;  *) TARGET_FILE="$DIR/$TARGET_FILE" ;;
    esac
  done
  DIR=$(dirname -- "$TARGET_FILE")
  (CDPATH= cd -- "$DIR" && pwd)
}
SCRIPT_DIR="$(get_script_dir)"

SITE="$SCRIPT_DIR/../share/pytopp/site"
export PYTHONPATH="${SITE}${PYTHONPATH:+:${PYTHONPATH}}"

TOOL="$SITE/pyopenms/pytopp/tools/pytopp_pyprophet_cli.py"

# Prefer uv (the CLI has: #!/usr/bin/env -S uv run --script)
if command -v uv >/dev/null 2>&1; then
  exec "$TOOL" "$@"
fi

# Fallback to python with a helpful message
if ! command -v python >/dev/null 2>&1; then
  cat <<'MSG'
ERROR: Neither `uv` nor `python` found on PATH.

Install uv (recommended) and re-run:
  macOS/Linux (curl):   curl -LsSf https://astral.sh/uv/install.sh | sh
  macOS (Homebrew):     brew install uv
  Linux (pipx):         pipx install uv

Or install Python and run:
  python "$TOOL" --help
MSG
  exit 127
fi

exec python "$TOOL" "$@"
]=])
  file(CHMOD "${_wrapper}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                GROUP_READ GROUP_EXECUTE
                WORLD_READ WORLD_EXECUTE)
  install(PROGRAMS "${_wrapper}" DESTINATION bin)
endif()

# --------------------------- Install payload -------------------------------
# 1) Install all package .py files (including sitecustomize.py; exclude CLI itself)
install(DIRECTORY "${_pytopp_site_dir}/"
        DESTINATION "share/pytopp/site"
        FILES_MATCHING
          PATTERN "*.py"
          PATTERN "__pycache__" EXCLUDE
          PATTERN "*/${_pytopp_cli_rel}" EXCLUDE)

# 2) Install the CLI with exec perms preserved
install(PROGRAMS "${_pytopp_cli_build}"
        DESTINATION "share/pytopp/site/pyopenms/pytopp/tools")
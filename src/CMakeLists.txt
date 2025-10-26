# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_src")

#------------------------------------------------------------------------------
# OpenSWATHAlgo lib
#------------------------------------------------------------------------------
add_subdirectory(openswathalgo)

#------------------------------------------------------------------------------
# OpenMS core lib
#------------------------------------------------------------------------------
add_subdirectory(openms)

#------------------------------------------------------------------------------
# pyOpenMS
#------------------------------------------------------------------------------
option(PYOPENMS "setup build system for pyOpenMS" OFF)
if(PYOPENMS)
  add_subdirectory(pyOpenMS)
endif()

if(WITH_GUI)
  #------------------------------------------------------------------------------
  # OpenMS gui lib
  #------------------------------------------------------------------------------
  add_subdirectory(openms_gui)
endif()

# --- pyTOPP: ship the whole package + an OS-appropriate wrapper ---------------
option(PYTOPP "Install pyTOPP CLI wrappers (Python tools via uv)" ON)

if(PYTOPP)
  # Source layout (your repo)
  # We copy ONLY the 'pytopp' subpackage to keep 'pyopenms' as a namespace.
  set(_pytopp_pkg_src "${CMAKE_SOURCE_DIR}/src/pyTOPP/src/pyopenms/pytopp")

  # Where we stage things in the build/install trees
  set(_pytopp_site_dir "${CMAKE_BINARY_DIR}/share/pytopp/site")
  set(_pytopp_pkg_dst_root "${_pytopp_site_dir}/pyopenms")
  file(MAKE_DIRECTORY "${_pytopp_pkg_dst_root}")

  # Copy the whole package: share/pytopp/site/pyopenms/pytopp/...
  file(COPY "${_pytopp_pkg_src}" DESTINATION "${_pytopp_pkg_dst_root}")

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
    set(_wrapper "${CMAKE_BINARY_DIR}/bin/pyTOPP-PyProphet.cmd")
    file(WRITE "${_wrapper}" [=[
@echo off
REM Resolve wrapper dir
set "WRAPDIR=%~dp0"
if "%WRAPDIR:~-1%"=="\" set "WRAPDIR=%WRAPDIR:~0,-1%"

REM Vendored site root (..\\share\\pytopp\\site)
set "SITE=%WRAPDIR%\..\share\pytopp\site"
set "PYTHONPATH=%SITE%;%PYTHONPATH%"

REM CLI inside vendored package
set "TOOL=%SITE%\pyopenms\pytopp\tools\pytopp_pyprophet_cli.py"

REM Prefer uv (script's shebang uses `uv run --script`)
where uv >NUL 2>&1
if %ERRORLEVEL% EQU 0 (
  "%TOOL%" %*
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
    set(_wrapper "${CMAKE_BINARY_DIR}/bin/pyTOPP-PyProphet")
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
  # 1) Install all package .py files (except the CLI; we install it as PROGRAM)
  install(DIRECTORY "${_pytopp_site_dir}/"
          DESTINATION "share/pytopp/site"
          FILES_MATCHING
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE
            PATTERN "*/${_pytopp_cli_rel}" EXCLUDE)

  # 2) Install the CLI with exec perms preserved
  install(PROGRAMS "${_pytopp_cli_build}"
          DESTINATION "share/pytopp/site/pyopenms/pytopp/tools")
endif()

#------------------------------------------------------------------------------
# TOPP tools
#------------------------------------------------------------------------------
add_subdirectory(topp)

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------
add_subdirectory(tests)

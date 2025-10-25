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

# --- pyTOPP: ship the script and create an OS-appropriate wrapper -------------
option(PYTOPP "Install pyTOPP CLI wrappers (Python tools via uv)" ON)

if(PYTOPP)
  # Source tool (your uv-run script)
  set(_pytopp_tool_src
      "${CMAKE_SOURCE_DIR}/src/pyTOPP/src/pyopenms/pytopp/tools/pytopp_pyprophet_cli.py")

  # Where we ship it inside the build/install layout
  set(_pytopp_rel_tool "share/pytopp/tools/pytopp_pyprophet_cli.py")
  set(_pytopp_tool_dst "${CMAKE_BINARY_DIR}/${_pytopp_rel_tool}")

  # Ensure dir and copy tool
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/share/pytopp/tools")
  file(COPY "${_pytopp_tool_src}" DESTINATION "${CMAKE_BINARY_DIR}/share/pytopp/tools")

  # Executable bit for the shipped tool on Unix
  if(UNIX)
    file(CHMOD "${_pytopp_tool_dst}"
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                  GROUP_READ GROUP_EXECUTE
                  WORLD_READ WORLD_EXECUTE)
  endif()

  # Wrapper per host OS
  if(WIN32)
    # -------- Windows wrapper (.cmd) ------------------------------------------
    set(_win_wrapper "${CMAKE_BINARY_DIR}/bin/pyTOPP-PyProphet.cmd")
    file(WRITE "${_win_wrapper}" [=[
@echo off
REM Resolve wrapper dir
set "WRAPDIR=%~dp0"
if "%WRAPDIR:~-1%"=="\" set "WRAPDIR=%WRAPDIR:~0,-1%"

REM Tool path relative to bin\ -> ..\share\pytopp\tools\...
set "TOOL=%WRAPDIR%\..\share\pytopp\tools\pytopp_pyprophet_cli.py"

REM Friendly check for uv (script’s shebang uses `uv run`)
where uv >NUL 2>&1
if %ERRORLEVEL% NEQ 0 (
  echo ERROR: "uv" is not installed, but this tool uses uv to auto-manage Python ^& dependencies.
  echo.
  echo Install one of the following and re-run this command:
  echo   Windows (PowerShell):  iwr https://astral.sh/uv/install.ps1 -UseB ^| iex
  echo   Windows (winget):      winget install --id AstralSoftware.UV
  echo.
  echo Advanced/offline: run with an existing Python instead:
  echo   python "%TOOL%" --help
  exit /b 127
)

REM Hand off to the shipped script (its shebang invokes `uv run`)
"%TOOL%" %*
]=])

    # Install Windows wrapper
    install(PROGRAMS "${_win_wrapper}" DESTINATION bin)

  else()
    # -------- POSIX wrapper (macOS/Linux) -------------------------------------
    set(_posix_wrapper "${CMAKE_BINARY_DIR}/bin/pyTOPP-PyProphet")
    file(WRITE "${_posix_wrapper}" [=[
#!/usr/bin/env sh
# Resolve path to this wrapper (works with symlinks)
get_script_dir() {
  TARGET_FILE="$0"
  while [ -L "$TARGET_FILE" ]; do
    DIR=$(dirname -- "$TARGET_FILE")
    TARGET_FILE=$(readlink -- "$TARGET_FILE")
    case "$TARGET_FILE" in
      /*) ;; # absolute
      *) TARGET_FILE="$DIR/$TARGET_FILE" ;;
    esac
  done
  DIR=$(dirname -- "$TARGET_FILE")
  (CDPATH= cd -- "$DIR" && pwd)
}
SCRIPT_DIR="$(get_script_dir)"
TOOL="$SCRIPT_DIR/../share/pytopp/tools/pytopp_pyprophet_cli.py"

# Friendly check for uv (script’s shebang uses `uv run`)
if ! command -v uv >/dev/null 2>&1; then
  cat <<'MSG'
ERROR: `uv` is not installed, but this tool uses uv to auto-manage Python & dependencies.

Install one of the following and re-run this command:
  macOS/Linux (curl):   curl -LsSf https://astral.sh/uv/install.sh | sh
  macOS (Homebrew):     brew install uv
  Linux (pipx):         pipx install uv

Advanced/offline: you can also run the script with your own Python:
  python "$TOOL" --help
MSG
  exit 127
fi

# Hand off to the shipped script (has: #!/usr/bin/env -S uv run --script)
exec "$TOOL" "$@"
]=])

    # Exec bit + install POSIX wrapper
    file(CHMOD "${_posix_wrapper}"
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                  GROUP_READ GROUP_EXECUTE
                  WORLD_READ WORLD_EXECUTE)
    install(PROGRAMS "${_posix_wrapper}" DESTINATION bin)
  endif()

  # Install the shipped tool itself
  install(PROGRAMS "${_pytopp_tool_dst}" DESTINATION share/pytopp/tools)
endif()


#------------------------------------------------------------------------------
# TOPP tools
#------------------------------------------------------------------------------
add_subdirectory(topp)

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------
add_subdirectory(tests)

# Copyright (c) 2002-present, OpenMS Inc.
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------
# 3.15 to enable new find policies for FindPython
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("pyTOPP" LANGUAGES C CXX)

# Put generated pyTOPP wrappers next to TOPP binaries in the top-level bin/
set(_PYTOPP_WRAPPER_BIN "${CMAKE_BINARY_DIR}/bin")
file(MAKE_DIRECTORY "${_PYTOPP_WRAPPER_BIN}")

# ─────────────────────────────────────────────────────────────────────────────
# Package layout
# We copy ONLY the 'pytopp' subpackage so that 'pyopenms' stays a namespace.
# Runtime bootstrapping happens via sitecustomize.py we generate below.
# ─────────────────────────────────────────────────────────────────────────────
set(_pytopp_pkg_src      "${CMAKE_CURRENT_SOURCE_DIR}/src/pyopenms/pytopp")
set(_pytopp_site_dir     "${CMAKE_CURRENT_BINARY_DIR}/share/pytopp/site")
set(_pytopp_pkg_dst_root "${_pytopp_site_dir}/pyopenms")
file(MAKE_DIRECTORY "${_pytopp_pkg_dst_root}")

# Copy the whole package: share/pytopp/site/pyopenms/pytopp/...
file(COPY "${_pytopp_pkg_src}" DESTINATION "${_pytopp_pkg_dst_root}")

# ─────────────────────────────────────────────────────────────────────────────
# sitecustomize.py
# Ensures our vendored 'pyopenms/pytopp' is discoverable even if a system
# pyopenms is installed. Python auto-loads this module when its directory
# is on sys.path.
# ─────────────────────────────────────────────────────────────────────────────
file(WRITE "${_pytopp_site_dir}/sitecustomize.py" [=[
# Auto-imported by Python's 'site' module if this directory is on sys.path.
# Extends/creates the 'pyopenms' namespace to include the vendored 'pytopp'.

import importlib
import importlib.machinery
import os
import sys
import types

def _ensure_pyopenms_namespace():
    site = os.path.dirname(os.path.abspath(__file__))
    pyo_dir = os.path.join(site, "pyopenms")

    try:
        m = importlib.import_module("pyopenms")
    except Exception:
        m = None

    if m is not None:
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

    m = types.ModuleType("pyopenms")
    m.__package__ = "pyopenms"
    m.__path__ = [pyo_dir]
    m.__spec__ = importlib.machinery.ModuleSpec("pyopenms", loader=None, is_package=True)
    sys.modules["pyopenms"] = m

_ensure_pyopenms_namespace()
]=])

# ─────────────────────────────────────────────────────────────────────────────
# Helper: add_pytopp_tool
#
# Usage (one line per tool):
#   add_pytopp_tool(
#       NAME         PyProphet
#       CLI_REL      pyopenms/pytopp/tools/pytopp_pyprophet_cli.py
#   )
#
# What this does:
#   • Marks the CLI script in the build-tree as executable on UNIX
#   • Generates a wrapper in bin/:   bin/<NAME>
#   • Wrapper sets PYTHONPATH to include our vendored site tree
#   • Prefers running via `uv run --script`, falls back to `python`
#   • Installs the CLI script (with exec perms preserved) under share/pytopp/site
#   • Installs the wrapper into bin/
#
# Developer notes:
#   • To add a new tool, drop a new CLI under:
#       src/pyTOPP/src/pyopenms/pytopp/tools/<your_cli>.py
#     Give it a shebang like: `#!/usr/bin/env -S uv run --script`
#   • Then add ONE call to add_pytopp_tool(NAME <NiceName> CLI_REL <relative_path>)
#     below in the "Register tools" section.
#   • No other changes are needed.
# ─────────────────────────────────────────────────────────────────────────────
function(add_pytopp_tool)
  set(options)
  set(oneValueArgs NAME CLI_REL)
  set(multiValueArgs)
  cmake_parse_arguments(PYT "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT PYT_NAME)
    message(FATAL_ERROR "add_pytopp_tool(): NAME is required")
  endif()
  if(NOT PYT_CLI_REL)
    message(FATAL_ERROR "add_pytopp_tool(): CLI_REL is required")
  endif()

  # Paths to the CLI inside the staged vendored site
  set(_cli_rel   "${PYT_CLI_REL}")
  set(_cli_build "${_pytopp_site_dir}/${_cli_rel}")

  # Make the CLI executable on Unix (so a `uv` shebang works)
  if(UNIX)
    file(CHMOD "${_cli_build}"
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                  GROUP_READ  GROUP_EXECUTE
                  WORLD_READ  WORLD_EXECUTE)
  endif()

  # Wrapper filename in bin/
  if(WIN32)
    set(_wrapper "${_PYTOPP_WRAPPER_BIN}/${PYT_NAME}.cmd")
    file(WRITE "${_wrapper}" "@echo off\n")
    file(APPEND "${_wrapper}" "setlocal\n\n")
    file(APPEND "${_wrapper}" "REM Resolve wrapper dir\n")
    file(APPEND "${_wrapper}" "set \"WRAPDIR=%~dp0\"\n")
    file(APPEND "${_wrapper}" "if \"%WRAPDIR:~-1%\"==\"\\\" set \"WRAPDIR=%WRAPDIR:~0,-1%\"\n\n")
    file(APPEND "${_wrapper}" "REM Vendored site root (..\\share\\pytopp\\site)\n")
    file(APPEND "${_wrapper}" "set \"SITE=%WRAPDIR%\\..\\share\\pytopp\\site\"\n")
    file(APPEND "${_wrapper}" "set \"PYTHONPATH=%SITE%;%PYTHONPATH%\"\n\n")
    file(APPEND "${_wrapper}" "REM CLI inside vendored package\n")
    file(APPEND "${_wrapper}" "set \"TOOL=%SITE%\\${_cli_rel}\"\n\n")
    file(APPEND "${_wrapper}" "REM Prefer uv explicitly (Windows doesn't honor shebang)\n")
    file(APPEND "${_wrapper}" "where uv >NUL 2>&1\n")
    file(APPEND "${_wrapper}" "if %ERRORLEVEL% EQU 0 (\n")
    file(APPEND "${_wrapper}" "  uv run --script \"%TOOL%\" %*\n")
    file(APPEND "${_wrapper}" "  exit /b %ERRORLEVEL%\n")
    file(APPEND "${_wrapper}" ")\n\n")
    file(APPEND "${_wrapper}" "REM Fallback to python with a helpful message\n")
    file(APPEND "${_wrapper}" "where python >NUL 2>&1\n")
    file(APPEND "${_wrapper}" "if %ERRORLEVEL% NEQ 0 (\n")
    file(APPEND "${_wrapper}" "  echo ERROR: Neither \"uv\" nor \"python\" found on PATH.\n")
    file(APPEND "${_wrapper}" "  echo.\n")
    file(APPEND "${_wrapper}" "  echo Install uv (recommended) and re-run:\n")
    file(APPEND "${_wrapper}" "  echo   PowerShell:  iwr https://astral.sh/uv/install.ps1 -UseB ^| iex\n")
    file(APPEND "${_wrapper}" "  echo   winget:      winget install --id AstralSoftware.UV\n")
    file(APPEND "${_wrapper}" "  echo.\n")
    file(APPEND "${_wrapper}" "  echo Or install Python and run:\n")
    file(APPEND "${_wrapper}" "  echo   python \"%TOOL%\" --help\n")
    file(APPEND "${_wrapper}" "  exit /b 127\n")
    file(APPEND "${_wrapper}" ")\n\n")
    file(APPEND "${_wrapper}" "python \"%TOOL%\" %*\n")

    # Install the wrapper
    install(PROGRAMS "${_wrapper}" DESTINATION bin)
  else()
    set(_wrapper "${_PYTOPP_WRAPPER_BIN}/${PYT_NAME}")
    file(WRITE "${_wrapper}" "#!/usr/bin/env sh\n")
    file(APPEND "${_wrapper}" "set -eu\n")
    file(APPEND "${_wrapper}" "get_script_dir() {\n")
    file(APPEND "${_wrapper}" "  TARGET_FILE=\"$0\"\n")
    file(APPEND "${_wrapper}" "  while [ -L \"$TARGET_FILE\" ]; do\n")
    file(APPEND "${_wrapper}" "    DIR=$(dirname -- \"$TARGET_FILE\")\n")
    file(APPEND "${_wrapper}" "    TARGET_FILE=$(readlink -- \"$TARGET_FILE\")\n")
    file(APPEND "${_wrapper}" "    case \"$TARGET_FILE\" in\n")
    file(APPEND "${_wrapper}" "      /*) ;;  *) TARGET_FILE=\"$DIR/$TARGET_FILE\" ;;\n")
    file(APPEND "${_wrapper}" "    esac\n")
    file(APPEND "${_wrapper}" "  done\n")
    file(APPEND "${_wrapper}" "  DIR=$(dirname -- \"$TARGET_FILE\")\n")
    file(APPEND "${_wrapper}" "  (CDPATH= cd -- \"$DIR\" && pwd)\n")
    file(APPEND "${_wrapper}" "}\n")
    file(APPEND "${_wrapper}" "SCRIPT_DIR=\"$(get_script_dir)\"\n\n")
    file(APPEND "${_wrapper}" "SITE=\"$SCRIPT_DIR/../share/pytopp/site\"\n")
    file(APPEND "${_wrapper}" "export PYTHONPATH=\"${SITE}\${PYTHONPATH:+:${PYTHONPATH}}\"\n\n")
    file(APPEND "${_wrapper}" "TOOL=\"$SITE/${_cli_rel}\"\n\n")
    file(APPEND "${_wrapper}" "# Prefer uv (the CLI may have: #!/usr/bin/env -S uv run --script)\n")
    file(APPEND "${_wrapper}" "if command -v uv >/dev/null 2>&1; then\n")
    file(APPEND "${_wrapper}" "  exec \"$TOOL\" \"$@\"\n")
    file(APPEND "${_wrapper}" "fi\n\n")
    file(APPEND "${_wrapper}" "# Fallback to python with a helpful message\n")
    file(APPEND "${_wrapper}" "if ! command -v python >/dev/null 2>&1; then\n")
    file(APPEND "${_wrapper}" "  cat <<'MSG'\n")
    file(APPEND "${_wrapper}" "ERROR: Neither \`uv\` nor \`python\` found on PATH.\n\n")
    file(APPEND "${_wrapper}" "Install uv (recommended) and re-run:\n")
    file(APPEND "${_wrapper}" "  macOS/Linux (curl):   curl -LsSf https://astral.sh/uv/install.sh | sh\n")
    file(APPEND "${_wrapper}" "  macOS (Homebrew):     brew install uv\n")
    file(APPEND "${_wrapper}" "  Linux (pipx):         pipx install uv\n\n")
    file(APPEND "${_wrapper}" "Or install Python and run:\n")
    file(APPEND "${_wrapper}" "  python \"$TOOL\" --help\n")
    file(APPEND "${_wrapper}" "MSG\n")
    file(APPEND "${_wrapper}" "  exit 127\n")
    file(APPEND "${_wrapper}" "fi\n\n")
    file(APPEND "${_wrapper}" "exec python \"$TOOL\" \"$@\"\n")

    file(CHMOD "${_wrapper}"
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                  GROUP_READ  GROUP_EXECUTE
                  WORLD_READ  WORLD_EXECUTE)

    # Install the wrapper
    install(PROGRAMS "${_wrapper}" DESTINATION bin)
  endif()

  # Install the CLI itself (preserve exec perms on UNIX)
  install(PROGRAMS "${_cli_build}"
          DESTINATION "share/pytopp/site/pyopenms/pytopp/tools")
endfunction()

# ─────────────────────────────────────────────────────────────────────────────
# Register tools here
#
# HOW TO ADD A NEW pyTOPP TOOL:
#   1) Create your CLI script under:
#        src/pyTOPP/src/pyopenms/pytopp/tools/<your_cli>.py
#      Suggested name pattern: "pytopp_<name>_cli.py"
#   2) Add ONE line below, e.g.:
#        add_pytopp_tool(NAME MyNewTool
#                        CLI_REL pyopenms/pytopp/tools/pytopp_mynewtool_cli.py)
#   3) Reconfigure & build. You'll get a wrapper:
#        build/.../bin/MyNewTool
#      and the CLI is installed under:
#        share/pytopp/site/pyopenms/pytopp/tools/<your_cli>.py
#   4) Dev tip: you can still run your CLI directly during development if you
#      export PYTHONPATH=src/pyTOPP/src (or add the tiny bootstrap shim).
# ─────────────────────────────────────────────────────────────────────────────

# Example: existing PyProphet tool
add_pytopp_tool(
  NAME    OpenSwathProphet
  CLI_REL pyopenms/pytopp/tools/pytopp_pyprophet_cli.py
)

# ─────────────────────────────────────────────────────────────────────────────
# Install the rest of the package payload (excluding tool scripts, which we
# install as PROGRAMS above to preserve executable bits).
# ─────────────────────────────────────────────────────────────────────────────
install(DIRECTORY "${_pytopp_site_dir}/"
        DESTINATION "share/pytopp/site"
        FILES_MATCHING
          PATTERN "*.py"
          PATTERN "__pycache__" EXCLUDE
          # Exclude ALL tool .py files; they are installed via add_pytopp_tool()
          PATTERN "pyopenms/pytopp/tools/*.py" EXCLUDE)


# ─────────────────────────────────────────────────────────────────────────────
# Testin pyTOPP
# ─────────────────────────────────────────────────────────────────────────────

# Toggle tests for this subproject
option(PYTOPP_ENABLE_TESTS "Enable pytest-based tests for pyTOPP" ON)

if(PYTOPP_ENABLE_TESTS)
  enable_testing()

  # Python interpreter (mirror pyOpenMS style variable if needed)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  if(NOT Python_EXECUTABLE)
    set(Python_EXECUTABLE ${Python3_EXECUTABLE})
  endif()

  # Where tests live and where they should run from
  set(_PYTOPP_TEST_WORKDIR "${CMAKE_CURRENT_SOURCE_DIR}")

  # Environment so tests can import:
  #  - source tree (live dev code)
  #  - build vendored site (sitecustomize.py + staged package)
  set(_PYTOPP_TEST_ENV
      "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/src;${CMAKE_CURRENT_BINARY_DIR}/share/pytopp/site;$ENV{PYTHONPATH}"
  )
  if(NOT WIN32)
    list(APPEND _PYTOPP_TEST_ENV "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
  endif()

  # Run all pytest tests in one go (fast feedback, but less granular on CDash)
  add_test(
    NAME pytopp_pytest_all
    COMMAND ${Python_EXECUTABLE} -m pytest -q tests
    WORKING_DIRECTORY ${_PYTOPP_TEST_WORKDIR}
  )
  set_tests_properties(pytopp_pytest_all PROPERTIES ENVIRONMENT "${_PYTOPP_TEST_ENV}")

  # -----------------------------------------------------------------------
  # Please add your unit test filenames here (under src/pyTOPP/tests/)
  # -----------------------------------------------------------------------
  set(pytopp_unittest_testfiles
    test_util.py
    test_ctd_support.py
  )

  # -----------------------------------------------------------------------
  # Please add your integration test filenames here (under src/pyTOPP/tests/)
  # -----------------------------------------------------------------------
  set(pytopp_integrationtest_testfiles
    test_integration_pytopp_pyprophet.py
  )

  # Create one CTest entry per unit test file
  foreach(t ${pytopp_unittest_testfiles})
    add_test(
      NAME "pytopp_unittest_${t}"
      COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/tests/${t}
      WORKING_DIRECTORY ${_PYTOPP_TEST_WORKDIR}
    )
    set_tests_properties("pytopp_unittest_${t}" PROPERTIES ENVIRONMENT "${_PYTOPP_TEST_ENV}")
  endforeach()

  # Create one CTest entry per integration test file
  foreach(t ${pytopp_integrationtest_testfiles})
    add_test(
      NAME "pytopp_integrationtest_${t}"
      COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/tests/${t}
      WORKING_DIRECTORY ${_PYTOPP_TEST_WORKDIR}
    )
    set_tests_properties("pytopp_integrationtest_${t}" PROPERTIES ENVIRONMENT "${_PYTOPP_TEST_ENV}")
  endforeach()

  # Optional: smoke the built wrapper in bin/ (prints --help). Safe to remove if you prefer.
  if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/bin/OpenSwathProphet")
    add_test(
      NAME pytopp_wrapper_OpenSwathProphet_help
      COMMAND "${CMAKE_CURRENT_BINARY_DIR}/bin/OpenSwathProphet" --help
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
    set_tests_properties(pytopp_wrapper_OpenSwathProphet_help
      PROPERTIES ENVIRONMENT "${_PYTOPP_TEST_ENV}")
  endif()

endif()
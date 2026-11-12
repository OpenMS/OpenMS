# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


project("OpenMS_src")

#------------------------------------------------------------------------------
# The libraries. A build against an installed OpenMS (OPENMS_USE_INSTALLED_LIBRARY)
# takes OpenSwathAlgo, OpenMS and OpenMS_CLI from that installation instead
# (cmake/find_installed_openms.cmake) and builds only the TOPP tools on top.
#------------------------------------------------------------------------------
if(NOT OPENMS_USE_INSTALLED_LIBRARY)
  #------------------------------------------------------------------------------
  # OpenSWATHAlgo lib
  #------------------------------------------------------------------------------
  add_subdirectory(openswathalgo)

  #------------------------------------------------------------------------------
  # OpenMS core lib
  #------------------------------------------------------------------------------
  add_subdirectory(openms)
  add_subdirectory(openms_cli)
endif()

#------------------------------------------------------------------------------
# Class-test framework lib (ClassTest macros, FuzzyStringComparator);
# standard-library-only; used by the class tests of all OpenMS libraries and
# by the FuzzyDiff tool
#------------------------------------------------------------------------------
add_subdirectory(testframework)

#------------------------------------------------------------------------------
# pyOpenMS
#------------------------------------------------------------------------------
# Note: pyOpenMS is optional and can be disabled with -DPYOPENMS=OFF
option(PYOPENMS "Build pyOpenMS Python bindings" OFF)
option(WITH_UV "Install Python dependencies with uv" ON)

# Optional: Python version to use when creating the pyOpenMS venv with uv
# Example: "3.11". Only used when WITH_UV=ON.
set(PYOPENMS_UV_PYTHON_VERSION "" CACHE STRING "Python version (e.g. 3.11) for creating pyOpenMS .venv when WITH_UV=ON")

if(PYOPENMS AND OPENMS_USE_INSTALLED_LIBRARY)
  message(FATAL_ERROR "PYOPENMS is not supported together with OPENMS_USE_INSTALLED_LIBRARY.")
endif()
if(PYOPENMS)
  if(WITH_UV)
    # Discover uv executable
    find_program(UV_EXECUTABLE uv)
    if(NOT UV_EXECUTABLE)
      message(FATAL_ERROR "uv executable not found. Install uv (https://docs.astral.sh/uv/) and ensure it's on PATH.")
    endif()
    message(STATUS "Found uv: ${UV_EXECUTABLE}")

    set(Python_FIND_VIRTUALENV FIRST)
    # We found uv; pyOpenMS will decide later whether to run it and where to place the .venv
  endif()

  message(STATUS "pyOpenMS Python environment: ${Python_EXECUTABLE}")

  # Build pyOpenMS Python extension modules directly as part of the main build
  # This ensures all CMake variables are naturally inherited
  add_subdirectory(pyOpenMS)

  # Optional: Package pyOpenMS wheel from already-built modules
  # Usage: cmake --build . --target pyopenms_wheel
  add_custom_target(pyopenms_wheel
    COMMAND ${Python_EXECUTABLE} -m build --wheel
      --skip-dependency-check --no-isolation
      --outdir ${CMAKE_BINARY_DIR}/pyopenms_wheels
      --config-setting=override=cmake.0.build_path=${CMAKE_BINARY_DIR}/py-build-cmake_cache/pyopenms
      --config-setting=override=cmake.0.options.PYOPENMS_USE_PREBUILT=ON
      --config-setting=override=cmake.0.options.PYOPENMS_PREBUILT_DIR=${CMAKE_BINARY_DIR}/pyOpenMS
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/pyOpenMS
    DEPENDS pyopenms
    COMMENT "Packaging pyOpenMS wheel from existing build"
  )
endif()

#------------------------------------------------------------------------------
# JavaScript bindings of the WebAssembly build (cmake/emscripten.cmake)
#------------------------------------------------------------------------------
if(EMSCRIPTEN AND OPENMS_WASM_BINDINGS AND NOT OPENMS_USE_INSTALLED_LIBRARY)
  add_subdirectory(openms_wasm)
endif()

#------------------------------------------------------------------------------
# TOPP tools
#------------------------------------------------------------------------------
# Declarations accumulate in the cache, so drop what an earlier configure of this build
# directory left behind before any tool is declared.
openms_reset_topp_tools()
if(BUILD_TOPP_TOOLS)
  add_subdirectory(topp)
endif()

#------------------------------------------------------------------------------
# The tool registry, generated from the openms_topp_tool() declarations of src/topp. It
# therefore describes exactly the tools this build produces: a build option that does not
# build a tool does not declare it either, so nothing has to filter the registry afterwards,
# at run time or at install time.
#------------------------------------------------------------------------------
#
# It is installed into TOOLS/ rather than beside the other share data so that the directory
# holds registry files and nothing else: ToolHandler reads every *.tsv it finds there.
#
# A build of the whole of OpenMS writes one registry, OpenMS.tsv, which belongs to the shared
# data. A build against an installed OpenMS library (OPENMS_USE_INSTALLED_LIBRARY) is one
# layer of a split installation: it writes the registry of the tools it builds,
# OpenMS-TOPP.tsv, installed with those tools (install component Applications), so each
# layer contributes its own file to a TOOLS directory they share (see
# OPENMS_INSTALL_DATA_PATH in src/openms/configh.cmake). The GUI layer (OpenMS-TOPP-GUI)
# adds OpenMS-GUI.tsv for its TOPP tools ExecutePipeline and ImageCreator.
#------------------------------------------------------------------------------
function(_openms_install_tool_registry file component)
  if(EXISTS "${file}")
    install(FILES "${file}"
      DESTINATION ${INSTALL_SHARE_DIR}/TOOLS
      COMPONENT ${component}
      PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ
    )
  endif()
endfunction()

if(NOT OPENMS_USE_INSTALLED_LIBRARY)
  openms_write_tool_registry("${OPENMS_TOOL_REGISTRY_BUILD_FILE}")
  _openms_install_tool_registry("${OPENMS_TOOL_REGISTRY_BUILD_FILE}" share)
  ## the registry file of a layered build of an earlier configure of this build directory
  file(REMOVE "${OPENMS_TOOL_REGISTRY_BUILD_DIR}/OpenMS-TOPP.tsv")
else()
  set(_registry_file "${OPENMS_TOOL_REGISTRY_BUILD_DIR}/OpenMS-TOPP.tsv")
  openms_write_tool_registry("${_registry_file}")
  _openms_install_tool_registry("${_registry_file}" Applications)
  ## the registry file of a full build of an earlier configure of this build directory
  file(REMOVE "${OPENMS_TOOL_REGISTRY_BUILD_FILE}")
endif()
## the GUI tools' registry that a layered build of this build directory wrote before the
## GUI moved to OpenMS-TOPP-GUI: it would list tools this build no longer has
file(REMOVE "${OPENMS_TOOL_REGISTRY_BUILD_DIR}/OpenMS-GUI.tsv")

#------------------------------------------------------------------------------
# Tests (they need the library targets of this build)
#------------------------------------------------------------------------------
if(NOT OPENMS_USE_INSTALLED_LIBRARY)
  add_subdirectory(tests)
endif()

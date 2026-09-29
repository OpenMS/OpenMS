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
# (cmake/find_installed_openms.cmake) and builds only the applications on top.
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
# TOPP tools (must come before openms_gui so TOPP_TOOLS cache is initialized
# before the GUI appends its tools to it)
#------------------------------------------------------------------------------
# Declarations accumulate in the cache, so drop what an earlier configure of this build
# directory left behind before any tool is declared.
openms_reset_topp_tools()
if(BUILD_TOPP_TOOLS)
  add_subdirectory(topp)
endif()
## the tools declared by src/topp, for the per-layer registry below
set(_openms_topp_layer_tools ${TOPP_TOOLS})

if(WITH_GUI)
  #------------------------------------------------------------------------------
  # OpenMS gui lib
  #------------------------------------------------------------------------------
  add_subdirectory(openms_gui)
else()
  # doc/doxygen/public references the GUI tools' documentation pages (the TOPP tool list, the
  # TOPPAS tutorial) whether or not the GUI is built, and what those pages describe does not
  # depend on building it. The part that does -- the command line help and INI documentation
  # TOPPDocumenter captures by running each tool -- is marked @if WITH_GUI in them (see
  # doc/CMakeLists.txt).
  openms_doc_path("${PROJECT_SOURCE_DIR}/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS")
endif()

#------------------------------------------------------------------------------
# The tool registry, generated from the openms_topp_tool() declarations of the two
# subdirectories above. It therefore describes exactly the tools this build produces: a
# build option that does not build a tool does not declare it either, so nothing has to
# filter the registry afterwards, at run time or at install time.
#------------------------------------------------------------------------------
#
# It is installed into TOOLS/ rather than beside the other share data so that the directory
# holds registry files and nothing else: ToolHandler reads every *.tsv it finds there.
#
# A build of the whole of OpenMS writes one registry, OpenMS.tsv, which belongs to the shared
# data. A build against an installed OpenMS library (OPENMS_USE_INSTALLED_LIBRARY) is one
# layer of a split installation: it writes one registry per group of tools it builds, installed
# with those tools, so each layer contributes its own file to a TOOLS directory they share
# (see OPENMS_INSTALL_DATA_PATH in src/openms/configh.cmake):
#   OpenMS-TOPP.tsv  the tools of src/topp          (install component Applications)
#   OpenMS-GUI.tsv   the TOPP tools of src/openms_gui (install component GUIApplications)
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
  ## the registry files of a layered build of an earlier configure of this build directory
  file(REMOVE "${OPENMS_TOOL_REGISTRY_BUILD_DIR}/OpenMS-TOPP.tsv" "${OPENMS_TOOL_REGISTRY_BUILD_DIR}/OpenMS-GUI.tsv")
else()
  set(_openms_gui_layer_tools ${TOPP_TOOLS})
  list(REMOVE_ITEM _openms_gui_layer_tools ${_openms_topp_layer_tools})
  foreach(_layer IN ITEMS TOPP GUI)
    string(TOLOWER "${_layer}" _layer_lower)
    set(_registry_file "${OPENMS_TOOL_REGISTRY_BUILD_DIR}/OpenMS-${_layer}.tsv")
    openms_write_tool_registry("${_registry_file}" TOOLS ${_openms_${_layer_lower}_layer_tools})
    if(_layer STREQUAL "TOPP")
      _openms_install_tool_registry("${_registry_file}" Applications)
    else()
      _openms_install_tool_registry("${_registry_file}" ${OPENMS_GUI_APPLICATIONS_COMPONENT})
    endif()
  endforeach()
  ## the registry file of a full build of an earlier configure of this build directory
  file(REMOVE "${OPENMS_TOOL_REGISTRY_BUILD_FILE}")
endif()

#------------------------------------------------------------------------------
# Tests (they need the library targets of this build)
#------------------------------------------------------------------------------
if(NOT OPENMS_USE_INSTALLED_LIBRARY)
  add_subdirectory(tests)
endif()

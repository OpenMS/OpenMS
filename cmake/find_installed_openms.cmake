# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Build the TOPP tools against an installed OpenMS library
# (OPENMS_USE_INSTALLED_LIBRARY=ON), used instead of cmake/cmake_findExternalLibs.cmake.
#
# The libraries OpenSwathAlgo, OpenMS and OpenMS_CLI, their headers and the shared
# data come from an installation of this same source tree (built with e.g.
# BUILD_TOPP_TOOLS=OFF); this build adds the TOPP tools. That is how OpenMS is packaged
# in layers, e.g. for Homebrew:
#   libopenms    the libraries, headers, CMake package and share/OpenMS
#   openms       the TOPP tools        (install component Applications)
#   openms-gui   OpenMS_GUI + TOPPView, TOPPAS, ... built from the OpenMS-TOPP-GUI
#                repository against the same installation
#
# Only what the TOPP tools themselves need is searched here: the dependencies of
# the libraries' public link interface are found by OpenMSConfig.cmake.
#------------------------------------------------------------------------------

if(NOT BUILD_TOPP_TOOLS)
  message(FATAL_ERROR "OPENMS_USE_INSTALLED_LIBRARY builds only the TOPP tools, so it needs BUILD_TOPP_TOOLS=ON.")
endif()
if(NOT "${PACKAGE_TYPE}" STREQUAL "none")
  message(FATAL_ERROR "OPENMS_USE_INSTALLED_LIBRARY does not support packaging (PACKAGE_TYPE=${PACKAGE_TYPE}).")
endif()

## The TOPP tools are compiled from the same sources as the library they link, so
## the installation has to be of this very version.
find_package(OpenMS ${OPENMS_PACKAGE_VERSION} EXACT CONFIG REQUIRED COMPONENTS CLI)
message(STATUS "Building against the installed OpenMS ${OpenMS_VERSION} in ${OpenMS_DIR}")
message(STATUS "  shared data of that installation: ${OPENMS_DATA_DIR}")

## options of the installed library that decide which applications exist
set(WITH_WNETALIGN ${OpenMS_WITH_WNETALIGN})

#------------------------------------------------------------------------------
# nlohmann::json: DIAuditor and OpenNuXL use it directly. It is a private dependency of
# libOpenMS, so the installation does not provide it; take it from where the library
# build does (src/openms/extern/CMakeLists.txt): the vendored copy, or an external one
# with USE_EXTERNAL_JSON.
option(USE_EXTERNAL_JSON "Use an external nlohmann::json library" OFF)
if(USE_EXTERNAL_JSON)
  find_package(nlohmann_json REQUIRED GLOBAL)
  target_compile_definitions(nlohmann_json::nlohmann_json INTERFACE JSON_USE_IMPLICIT_CONVERSIONS=0)
else()
  ## header-only: the include directory of the vendored copy, as its own CMakeLists.txt
  ## (which expects the setup of src/openms/extern/CMakeLists.txt) exports it
  add_library(nlohmann_json INTERFACE)
  add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)
  target_include_directories(nlohmann_json SYSTEM INTERFACE
    "${OPENMS_HOST_DIRECTORY}/src/openms/extern/nlohmann_json/include")
endif()

#------------------------------------------------------------------------------
# Boost (header-only): a few TOPP tools (FileInfo, FileMerger, TextExporter, ...; the
# list in src/topp/CMakeLists.txt) use it directly. It is a private
# dependency of libOpenMS, so the installation does not provide it; find it as the
# library build does.
find_boost()


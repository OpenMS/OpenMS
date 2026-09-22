# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(OpenNuXL)

# allow building separately from OpenMS
if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED)
endif()

# The whole NuXL pipeline. It used to live in OpenMS/ANALYSIS/NUXL, but nothing outside
# this tool ever used it, so it is kept next to the tool instead of in the library:
#   include/  the NuXL headers        source/  their implementations
#   test/     their class tests, built by the class-test project (see test/CMakeLists.txt)
set(tool_headers
  include/NuXLAnnotateAndLocate.h
  include/NuXLAnnotatedHit.h
  include/NuXLConstants.h
  include/NuXLDeisotoper.h
  include/NuXLFDR.h
  include/NuXLFragmentAdductDefinition.h
  include/NuXLFragmentAnnotationHelper.h
  include/NuXLFragmentIonGenerator.h
  include/NuXLMarkerIonExtractor.h
  include/NuXLModificationsGenerator.h
  include/NuXLParameterParsing.h
  include/NuXLPresets.h
  include/NuXLReport.h
)
set(tool_sources
  source/NuXLAnnotateAndLocate.cpp
  source/NuXLDeisotoper.cpp
  source/NuXLFDR.cpp
  source/NuXLFragmentAdductDefinition.cpp
  source/NuXLFragmentAnnotationHelper.cpp
  source/NuXLFragmentIonGenerator.cpp
  source/NuXLMarkerIonExtractor.cpp
  source/NuXLModificationsGenerator.cpp
  source/NuXLParameterParsing.cpp
  source/NuXLPresets.cpp
  source/NuXLReport.cpp
)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources} ${tool_headers})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI)
# NuXLPresets.cpp parses the preset JSON. nlohmann::json is a PRIVATE dependency of
# libOpenMS, so it does not arrive through OpenMS_CLI and has to be linked here.
target_link_libraries(${PROJECT_NAME} nlohmann_json::nlohmann_json)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

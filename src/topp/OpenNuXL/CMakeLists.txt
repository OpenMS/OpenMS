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
# this tool ever used it, so it is kept next to the tool instead of in the library.
# The class tests in src/tests/class_tests/openms build the sources they cover straight
# from this folder; see the OpenNuXL block in that directory's CMakeLists.txt.
set(tool_sources
  NuXLAnnotateAndLocate.cpp
  NuXLAnnotateAndLocate.h
  NuXLAnnotatedHit.h
  NuXLConstants.h
  NuXLDeisotoper.cpp
  NuXLDeisotoper.h
  NuXLFDR.cpp
  NuXLFDR.h
  NuXLFragmentAdductDefinition.cpp
  NuXLFragmentAdductDefinition.h
  NuXLFragmentAnnotationHelper.cpp
  NuXLFragmentAnnotationHelper.h
  NuXLFragmentIonGenerator.cpp
  NuXLFragmentIonGenerator.h
  NuXLMarkerIonExtractor.cpp
  NuXLMarkerIonExtractor.h
  NuXLModificationsGenerator.cpp
  NuXLModificationsGenerator.h
  NuXLParameterParsing.cpp
  NuXLParameterParsing.h
  NuXLPresets.cpp
  NuXLPresets.h
  NuXLReport.cpp
  NuXLReport.h
)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources})
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI)
# NuXLPresets.cpp parses the preset JSON. nlohmann::json is a PRIVATE dependency of
# libOpenMS, so it does not arrive through OpenMS_CLI and has to be linked here.
target_link_libraries(${PROJECT_NAME} nlohmann_json::nlohmann_json)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

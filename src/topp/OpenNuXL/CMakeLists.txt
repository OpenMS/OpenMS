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
# nlohmann::json is a PRIVATE dependency of libOpenMS, so the OpenMS package does not
# carry it: a standalone build has to look it up itself (NuXLPresets.cpp needs it).
if(NOT TARGET nlohmann_json::nlohmann_json)
  find_package(nlohmann_json CONFIG REQUIRED)
endif()
# Boost is a PRIVATE dependency of libOpenMS for the same reason, so it does not arrive
# through OpenMS_CLI either. Only the header-only target is needed: OpenNuXL.cpp uses
# Boost.Range and NuXLReport.cpp Boost.Math, neither of which is a compiled component.
if(NOT TARGET Boost::boost)
  find_package(Boost 1.81 REQUIRED)
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
# NuXLPresets.cpp parses the preset JSON; nlohmann::json does not arrive through OpenMS_CLI.
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI nlohmann_json::nlohmann_json Boost::boost)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

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

# Search/post-processing steps that only OpenNuXL uses. Kept next to the tool rather than in
# the library: they encode this pipeline's scoring, annotation and reporting decisions.
# The reusable NuXL parts (NuXLFragmentAdductDefinition, NuXLMarkerIonExtractor,
# NuXLModificationsGenerator, NuXLParameterParsing, NuXLPresets, NuXLReport) stay in
# OpenMS/ANALYSIS/NUXL, where their class tests can reach them.
set(tool_sources
  NuXLAnnotateAndLocate.cpp
  NuXLAnnotateAndLocate.h
  NuXLAnnotatedHit.h
  NuXLConstants.h
  NuXLDeisotoper.cpp
  NuXLDeisotoper.h
  NuXLFDR.cpp
  NuXLFDR.h
  NuXLFragmentAnnotationHelper.cpp
  NuXLFragmentAnnotationHelper.h
  NuXLFragmentIonGenerator.cpp
  NuXLFragmentIonGenerator.h
)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources})
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

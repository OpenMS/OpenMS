# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(MS1LabeledWorkflow)

# allow building separately from OpenMS
if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED)
endif()

# The channel-ratio quantification of this workflow. Kept next to the tool rather than in the
# library: its aggregation rules are the ones this workflow defines.
set(tool_sources
  MS1LabeledFAIMS.cpp
  MS1LabeledFAIMS.h
  MS1LabeledRatioQuantifier.cpp
  MS1LabeledRatioQuantifier.h
)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources})
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

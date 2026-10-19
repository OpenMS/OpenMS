# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(MetaProSIP)

# allow building separately from OpenMS
if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED)
endif()

# The pseudo features of identifications and spectra without a feature, kept next to the tool
# (with a class test in test/).
set(tool_sources
  MetaProSIPPseudoFeatures.cpp
  MetaProSIPPseudoFeatures.h
)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources})
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

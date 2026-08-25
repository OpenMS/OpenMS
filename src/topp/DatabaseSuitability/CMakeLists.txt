# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Tom Waschischeck $
# $Authors: Tom Waschischeck $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(DatabaseSuitability)

# allow building separately from OpenMS
if(NOT TARGET OpenMS)
  find_package(OpenMS)
endif()

# Classes which are specific to this tool live next to it instead of in the
# OpenMS library. They are compiled into a small static helper library
# (<ToolName>_lib) so that their class tests in src/tests/class_tests/openms
# can link against them as well (see executables.cmake there).
add_library(${PROJECT_NAME}_lib STATIC
  DBSuitability.h
  DBSuitability.cpp
)
openms_add_library_compiler_flags(${PROJECT_NAME}_lib)
target_include_directories(${PROJECT_NAME}_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# DBSuitability uses boost::regex directly
target_link_libraries(${PROJECT_NAME}_lib PUBLIC OpenMS Boost::regex)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
openms_add_executable_compiler_flags(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_lib)

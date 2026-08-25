# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Sven Nahnsen $
# $Authors: Sven Nahnsen $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(DecoyDatabase)

# allow building separately from OpenMS
if(NOT TARGET OpenMS)
  find_package(OpenMS)
endif()

# Classes which are specific to this tool live next to it instead of in the
# OpenMS library. They are compiled into a small static helper library
# (<ToolName>_lib) so that the tool's class tests (tests/) can link against
# them as well.
add_library(${PROJECT_NAME}_lib STATIC
  NeighborSeq.h
  NeighborSeq.cpp
)
openms_add_library_compiler_flags(${PROJECT_NAME}_lib)
target_include_directories(${PROJECT_NAME}_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME}_lib PUBLIC OpenMS)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
openms_add_executable_compiler_flags(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_lib Boost::regex)

# Class tests for the tool-specific classes above live with the tool (tests/);
# see the openms_add_tool_class_test macro in ../CMakeLists.txt.
openms_add_tool_class_test(NeighborSeq_test ${PROJECT_NAME}_lib)

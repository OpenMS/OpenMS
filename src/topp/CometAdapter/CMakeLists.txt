# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Chris Bielow $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(CometAdapter)

# allow building separately from OpenMS
if(NOT TARGET OpenMS)
  find_package(OpenMS)
endif()

# Classes which are specific to this tool live next to it instead of in the
# OpenMS library. They are compiled into a small static helper library
# (<ToolName>_lib) so that the tool's class tests (tests/) can link against
# them as well.
add_library(${PROJECT_NAME}_lib STATIC
  CometModification.h
  CometModification.cpp
  CometNativeIDRemapper.h
)
openms_add_library_compiler_flags(${PROJECT_NAME}_lib)
target_include_directories(${PROJECT_NAME}_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME}_lib PUBLIC OpenMS)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
openms_add_executable_compiler_flags(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_lib)

#------------------------------------------------------------------------------
# Class tests for the tool-specific classes above. Since src/topp is built on
# top of libOpenMS, these tests live with the tool (in tests/) instead of in
# src/tests/class_tests/openms. They honor the same switches as the library
# class tests. Note: ENABLE_CLASS_TESTING is declared in src/tests, which is
# processed after src/topp - hence the 'NOT DEFINED' fallback to its default ON.
if((NOT DEFINED ENABLE_CLASS_TESTING OR ENABLE_CLASS_TESTING) AND "${PACKAGE_TYPE}" STREQUAL "none")
  foreach(_tool_test CometModification_test CometNativeIDRemapper_test)
    add_executable(${_tool_test} tests/${_tool_test}.cpp)
    openms_add_executable_compiler_flags(${_tool_test})
    target_link_libraries(${_tool_test} ${PROJECT_NAME}_lib)
    # keep the test binaries out of the bin/ directory holding the TOPP tools
    set_target_properties(${_tool_test} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    add_test(NAME ${_tool_test} COMMAND $<TARGET_FILE:${_tool_test}>)
    if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set_target_properties(${_tool_test} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
    endif()
  endforeach()
endif()

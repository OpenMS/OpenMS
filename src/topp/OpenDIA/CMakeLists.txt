# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(OpenDIA)

if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
endif()

# Compile the orchestration once for both the executable and private boundary test.
# This target and its header are deliberately neither installed nor exported.
add_library(OpenDIA_impl STATIC source/OpenDIA.cpp include/OpenDIA.h
  source/OpenDIARecords.cpp include/OpenDIARecords.h)
target_include_directories(OpenDIA_impl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(OpenDIA_impl PUBLIC OpenMS::OpenMS_CLI)
if(COMMAND openms_add_library_compiler_flags)
  openms_add_library_compiler_flags(OpenDIA_impl)
endif()

add_executable(OpenDIA OpenDIA.cpp)
target_link_libraries(OpenDIA PRIVATE OpenDIA_impl)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(OpenDIA)
endif()

# In-tree tests are discovered later by src/tests/class_tests/openms.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(CTest)
  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()
endif()

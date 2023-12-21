# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# This file is called in the following scenarios
# 1. scenario: this file is at the root of a `cmake ..` call.
#       In this case we assume this project is compiled for packaging or testing.
# 2. scenario: another project is adding this directory via "find_package (TDL REQUIRED PATHS lib/tool_description_lib)".
#       In this case we don't want any tests to be included.

cmake_minimum_required (VERSION 3.15)

# check if this is the root project (
if (NOT DEFINED PROJECT_NAME)
    set (CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use")
    set (CMAKE_CXX_STANDARD_REQUIRED ON)
    set (CMAKE_CXX_EXTENSIONS OFF)
    set (TDL_ROOT_PROJECT TRUE)
endif ()

include (cmake/version.cmake)

project (tdl
         LANGUAGES CXX
         VERSION "${TDL_VERSION}"
         DESCRIPTION "tdl -- tool description library")

include (tdl-config.cmake)

# Enable and include testing if this is a root directory
if (${TDL_ROOT_PROJECT})
    enable_testing ()
    file (GLOB TEST_CPP_FILES
          LIST_DIRECTORIES false
          RELATIVE ${PROJECT_SOURCE_DIR}
          CONFIGURE_DEPENDS src/test_tdl/*.cpp)
    add_executable (test_tdl ${TEST_CPP_FILES})
    target_link_libraries (test_tdl tdl::tdl)
    add_test (NAME test_tdl COMMAND test_tdl)
endif ()

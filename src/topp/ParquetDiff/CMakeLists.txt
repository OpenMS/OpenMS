# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(ParquetDiff)

if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
endif()

# ParquetTableComparator live in libOpenMS, not here, although only this tool uses them: they
# are built on Arrow/Parquet, which are private dependencies of OpenMS. Kept in the library
# behind OpenMS-typed interfaces, they leave this tool with nothing to link but OpenMS, in-tree
# and against an installed OpenMS alike.
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

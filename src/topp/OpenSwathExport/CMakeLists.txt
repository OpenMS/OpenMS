# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(OpenSwathExport)

if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
endif()

# Arrow and Parquet are private dependencies of OpenMS. Resolve them here too
# so that this directory can also be built against an installed OpenMS SDK.
if(NOT OPENMS_ARROW_TARGET OR NOT OPENMS_PARQUET_TARGET)
  find_package(Arrow CONFIG REQUIRED)
  find_package(Parquet CONFIG REQUIRED)
  if(ARROW_USE_STATIC AND TARGET Arrow::arrow_static)
    set(OPENMS_ARROW_TARGET Arrow::arrow_static)
  elseif(TARGET Arrow::arrow_shared)
    set(OPENMS_ARROW_TARGET Arrow::arrow_shared)
  else()
    set(OPENMS_ARROW_TARGET Arrow::arrow_static)
  endif()
  if(ARROW_USE_STATIC AND TARGET Parquet::parquet_static)
    set(OPENMS_PARQUET_TARGET Parquet::parquet_static)
  elseif(TARGET Parquet::parquet_shared)
    set(OPENMS_PARQUET_TARGET Parquet::parquet_shared)
  else()
    set(OPENMS_PARQUET_TARGET Parquet::parquet_static)
  endif()
endif()

set(tool_sources
  source/OpenSwathMatrixExporter.cpp
  source/OpenSwathParquetExporter.cpp
  source/OpenSwathResultsExporter.cpp
)
set(tool_headers
  include/OpenSwathMatrixExporter.h
  include/OpenSwathParquetExporter.h
  include/OpenSwathResultsExporter.h
)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources} ${tool_headers})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET})
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

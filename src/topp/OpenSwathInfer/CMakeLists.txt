# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(OpenSwathInfer)

if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
endif()

set(tool_sources
  source/OpenSwathGeneInference.cpp
  source/OpenSwathPeptideInference.cpp
  source/OpenSwathPeptidoformInference.cpp
  source/OpenSwathProteinInference.cpp
)
set(tool_headers
  include/OpenSwathGeneInference.h
  include/OpenSwathPeptideInference.h
  include/OpenSwathPeptidoformInference.h
  include/OpenSwathProteinInference.h
)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources} ${tool_headers})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI )
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

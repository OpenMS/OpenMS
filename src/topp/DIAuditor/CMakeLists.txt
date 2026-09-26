# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg, David L. Tabb $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(DIAuditor)

if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
endif()
# nlohmann::json is a PRIVATE dependency of libOpenMS, so the OpenMS package does not
# carry it: DIAQCMetrics.cpp writes mzQC with it and looks it up itself.
if(NOT TARGET nlohmann_json::nlohmann_json)
  find_package(nlohmann_json CONFIG REQUIRED)
endif()

# The DIA metrics are used by this tool only, so they are kept next to it:
#   include/ and source/  DIAQCMetrics        test/  its class test (see test/CMakeLists.txt)
set(tool_sources
  source/DIAQCMetrics.cpp
)
set(tool_headers
  include/DIAQCMetrics.h
)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources} ${tool_headers})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI nlohmann_json::nlohmann_json)
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

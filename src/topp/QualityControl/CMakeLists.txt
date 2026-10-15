# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
endif()
project(QualityControl)

if(NOT TARGET OpenMS::OpenMS_CLI)
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
endif()

set(tool_sources
  source/Contaminants.cpp
  source/FWHM.cpp
  source/FragmentMassError.cpp
  source/MQEvidenceExporter.cpp
  source/MQExporterHelper.cpp
  source/MQMsmsExporter.cpp
  source/Ms2SpectrumStats.cpp
  source/MzCalibration.cpp
  source/PeptideMass.cpp
  source/RTAlignment.cpp
)
set(tool_headers
  include/Contaminants.h
  include/FWHM.h
  include/FragmentMassError.h
  include/MQEvidenceExporter.h
  include/MQExporterHelper.h
  include/MQMsmsExporter.h
  include/Ms2SpectrumStats.h
  include/MzCalibration.h
  include/PeptideMass.h
  include/RTAlignment.h
)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources} ${tool_headers})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} OpenMS::OpenMS_CLI )
if(COMMAND openms_add_executable_compiler_flags)
  openms_add_executable_compiler_flags(${PROJECT_NAME})
endif()

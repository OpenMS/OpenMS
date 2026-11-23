# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(NOT COMMAND openms_add_tool_class_test)
  message(FATAL_ERROR
    "${CMAKE_CURRENT_SOURCE_DIR} is added by src/tests/class_tests/openms, which provides "
    "openms_add_tool_class_test(). It cannot be configured on its own.")
endif()

openms_add_tool_class_test(Contaminants_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/Contaminants.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

openms_add_tool_class_test(FWHM_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/FWHM.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

openms_add_tool_class_test(FragmentMassError_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/FragmentMassError.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

openms_add_tool_class_test(MQEvidenceExporter_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/MQEvidenceExporter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/MQExporterHelper.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

openms_add_tool_class_test(MQMsmsExporter_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/MQExporterHelper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/MQMsmsExporter.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

openms_add_tool_class_test(Ms2SpectrumStats_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/Ms2SpectrumStats.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

openms_add_tool_class_test(MzCalibration_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/MzCalibration.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

openms_add_tool_class_test(PeptideMass_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/PeptideMass.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

openms_add_tool_class_test(RTAlignment_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/RTAlignment.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

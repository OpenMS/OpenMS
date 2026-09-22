# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(NOT COMMAND openms_add_tool_class_test)
  message(FATAL_ERROR
    "${CMAKE_CURRENT_SOURCE_DIR} is added by src/tests/class_tests/openms, which provides "
    "# Match the pre-relocation class-test gate.
if(DISABLE_OPENSWATH)
  return()
endif()

openms_add_tool_class_test(). It cannot be configured on its own.")
endif()

openms_add_tool_class_test(OpenSwathExport_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathMatrixExporter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathParquetExporter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathResultsExporter.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(OpenSwathExport_test ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET})

openms_add_tool_class_test(OpenSwathMatrixExporter_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathMatrixExporter.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(OpenSwathMatrixExporter_test ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET})

openms_add_tool_class_test(OpenSwathParquetExporter_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathParquetExporter.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(OpenSwathParquetExporter_test ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET})

openms_add_tool_class_test(OpenSwathResultsExporter_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathResultsExporter.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(OpenSwathResultsExporter_test ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET})

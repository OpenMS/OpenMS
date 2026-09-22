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

openms_add_tool_class_test(XIPMParquetConsumer_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/XIPMParquetConsumer.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(XIPMParquetConsumer_test ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET})

openms_add_tool_class_test(XIPMParquetRoundTrip_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/XIPMParquetConsumer.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(XIPMParquetRoundTrip_test ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET})

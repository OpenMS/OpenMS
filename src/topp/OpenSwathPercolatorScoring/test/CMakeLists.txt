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

# Match the pre-relocation class-test gate.
if(DISABLE_OPENSWATH)
  return()
endif()

# Use the same SQLite provider as OpenMS for in-tree builds.
if(TARGET unofficial::sqlite3::sqlite3)
  set(tool_sqlite unofficial::sqlite3::sqlite3)
elseif(TARGET SQLite::SQLite3)
  set(tool_sqlite SQLite::SQLite3)
elseif(TARGET SQLiteCpp)
  set(tool_sqlite SQLiteCpp ${OPENMS_SQLITECPP_EXTRA_LIBS})
elseif(TARGET sqlite3)
  set(tool_sqlite sqlite3)
else()
  find_package(SQLite3 REQUIRED)
  set(tool_sqlite SQLite::SQLite3)
endif()

openms_add_tool_class_test(OpenSwathOSWParquetReader_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathOSWParquetReader.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(OpenSwathOSWParquetReader_test ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET} ${tool_sqlite})

openms_add_tool_class_test(OpenSwathPercolatorScoring_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathOSWParquetReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathPercolatorScoring.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(OpenSwathPercolatorScoring_test ${OPENMS_ARROW_TARGET} ${OPENMS_PARQUET_TARGET} ${tool_sqlite})

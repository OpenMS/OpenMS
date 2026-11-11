# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# $Maintainer: Karnajeet Gosavi $

cmake_minimum_required(VERSION 3.24)

if(NOT IS_ABSOLUTE "${BINARY_DIR}")
  message(FATAL_ERROR "BINARY_DIR must name an absolute test output directory.")
endif()

if(DEFINED CASE)
  # Model CMake's discovery outputs, without requiring each CMake/MSVC version.
  set(OPENMS_IMPORT_STD OFF)
  include("${CMAKE_CURRENT_LIST_DIR}/../../../cmake/import_std.cmake")
  set(OPENMS_IMPORT_STD ON)
  set(CMAKE_CXX_COMPILER_ID MSVC)
  if(CASE STREQUAL "disabled")
    set(OPENMS_IMPORT_STD OFF)
    set(CMAKE_CXX_COMPILER_ID GNU)
  elseif(CASE STREQUAL "legacy_supported")
    set(CMAKE_VERSION 3.31.10)
    set(CMAKE_CXX_COMPILER_IMPORT_STD 23)
  elseif(CASE STREQUAL "legacy_stale")
    set(CMAKE_VERSION 3.31.10)
    set(CMAKE_CXX_COMPILER_IMPORT_STD "")
    set(CMAKE_CXX23_COMPILER_IMPORT_STD_NOT_FOUND_MESSAGE
      "Experimental import std support not enabled when detecting toolchain")
  elseif(CASE STREQUAL "modern_supported")
    set(CMAKE_VERSION 4.5)
    unset(CMAKE_CXX_COMPILER_IMPORT_STD)
    set(CMAKE_CXX_STDLIB_MODULES_JSON "${BINARY_DIR}/modules.json")
  elseif(CASE STREQUAL "modern_missing")
    set(CMAKE_VERSION 4.5)
    unset(CMAKE_CXX_COMPILER_IMPORT_STD)
    unset(CMAKE_CXX_STDLIB_MODULES_JSON)
    set(CMAKE_CXX_COMPILER_IMPORT_STD_ERROR_MESSAGE "Could not find modules.json resource")
  else()
    message(FATAL_ERROR "Unknown test case: ${CASE}")
  endif()
  openms_validate_import_std()
  return()
endif()

file(MAKE_DIRECTORY "${BINARY_DIR}")
file(WRITE "${BINARY_DIR}/modules.json" [=[
{"version":1,"revision":1,"modules":[{"logical-name":"std","source-path":"std.ixx","is-std-library":true}]}
]=])
foreach(_case IN ITEMS disabled legacy_supported legacy_stale modern_supported modern_missing)
  execute_process(COMMAND "${CMAKE_COMMAND}" "-DCASE=${_case}" "-DBINARY_DIR=${BINARY_DIR}"
    -P "${CMAKE_CURRENT_LIST_FILE}"
    RESULT_VARIABLE _result OUTPUT_VARIABLE _stdout ERROR_VARIABLE _stderr)
  # CMake wraps fatal diagnostics to the terminal width.
  string(REGEX REPLACE "[ \t\r\n]+" " " _stderr "${_stderr}")
  if(_case STREQUAL "legacy_stale" OR _case STREQUAL "modern_missing")
    if(_result EQUAL 0 OR NOT _stderr MATCHES "fresh build directory")
      message(FATAL_ERROR "${_case}: expected a failure with fresh-build guidance.\n${_stdout}\n${_stderr}")
    endif()
    if(_case STREQUAL "legacy_stale")
      set(_reason "not enabled when detecting toolchain")
    else()
      set(_reason "Could not find modules.json resource")
    endif()
    if(NOT _stderr MATCHES "${_reason}")
      message(FATAL_ERROR "${_case}: compiler-detection reason was lost.\n${_stderr}")
    endif()
  elseif(NOT _result EQUAL 0)
    message(FATAL_ERROR "${_case}: unexpected failure.\n${_stdout}\n${_stderr}")
  endif()
endforeach()

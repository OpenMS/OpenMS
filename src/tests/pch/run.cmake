# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Timo Sachsenberg $

cmake_minimum_required(VERSION 3.24)

if(NOT BINARY_DIR)
  message(FATAL_ERROR "Set BINARY_DIR to a scratch directory")
endif()
get_filename_component(BINARY_DIR "${BINARY_DIR}" ABSOLUTE)
file(MAKE_DIRECTORY "${BINARY_DIR}")
if(NOT GENERATOR)
  set(GENERATOR Ninja)
endif()
set(_configure -G "${GENERATOR}" -DCMAKE_BUILD_TYPE=Release)
foreach(_variable IN ITEMS CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_MAKE_PROGRAM
    CMAKE_CXX_COMPILER_LAUNCHER)
  if(DEFINED ${_variable})
    list(APPEND _configure "-D${_variable}=${${_variable}}")
  endif()
endforeach()

function(run step)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE _result OUTPUT_VARIABLE _out ERROR_VARIABLE _err)
  file(WRITE "${BINARY_DIR}/${step}.log" "${_out}${_err}")
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "${step} failed (${_result}):\n${_out}${_err}")
  endif()
endfunction()

# Reconfigure the same tree to catch stale PCH requirements when opting out.
foreach(_mode IN ITEMS OFF ON OFF)
  run(configure-${_mode} "${CMAKE_COMMAND}" -S "${CMAKE_CURRENT_LIST_DIR}"
    -B "${BINARY_DIR}/build" ${_configure} "-DOPENMS_USE_PCH=${_mode}")
  run(build-${_mode} "${CMAKE_COMMAND}" --build "${BINARY_DIR}/build"
    --config Release --parallel 2 --verbose)
  run(test-${_mode} "${CMAKE_CTEST_COMMAND}" --test-dir "${BINARY_DIR}/build"
    -C Release --output-on-failure)
  run(install-${_mode} "${CMAKE_COMMAND}" --install "${BINARY_DIR}/build"
    --config Release --prefix "${BINARY_DIR}/install")
  file(READ "${BINARY_DIR}/install/lib/cmake/PCHProbe/PCHProbeTargets.cmake" _export)
  if(_export MATCHES "PRECOMPILE_HEADERS|cmake_pch")
    message(FATAL_ERROR "Installed target exports build-specific PCH settings")
  endif()
endforeach()

# Test the MSVC /Yc bypass with direct options and response files. /Yu and
# ordinary compilations must still reach the compiler through the cache.
set(_trace "${BINARY_DIR}/launcher-trace.txt")
set(_mock "${CMAKE_CURRENT_LIST_DIR}/launcher.cmake")
set(_cache "${CMAKE_COMMAND}" -DROLE=cache "-DTRACE=${_trace}" -P "${_mock}" --)
list(LENGTH _cache _cache_count)
set(_compiler "${CMAKE_COMMAND}" -DROLE=compiler "-DTRACE=${_trace}" -P "${_mock}" --)
file(WRITE "${BINARY_DIR}/create pch.rsp" "/Yc\"header with spaces.hxx\" /Fp\"local.pch\"")
foreach(_case IN ITEMS create response use plain)
  if(_case STREQUAL "create")
    set(_args "/Ycheader with spaces.hxx")
    set(_expected "compiler\n[/Ycheader with spaces.hxx]\n")
  elseif(_case STREQUAL "response")
    set(_args "@${BINARY_DIR}/create pch.rsp")
    set(_expected "compiler\n[@${BINARY_DIR}/create pch.rsp]\n")
  elseif(_case STREQUAL "use")
    set(_args "/Yuheader with spaces.hxx")
    set(_expected "cache\ncompiler\n[/Yuheader with spaces.hxx]\n")
  else()
    set(_args "path with spaces.cpp" "-DVALUES=a\\;b")
    set(_expected "cache\ncompiler\n[path with spaces.cpp]\n[-DVALUES=a;b]\n")
  endif()
  file(WRITE "${_trace}" "")
  execute_process(COMMAND "${CMAKE_COMMAND}" "-DOPENMS_LAUNCHER_COUNT=${_cache_count}"
    -P "${CMAKE_CURRENT_LIST_DIR}/../../../cmake/msvc_pch_launcher.cmake" --
    ${_cache} ${_compiler} ${_args}
    RESULT_VARIABLE _result OUTPUT_VARIABLE _out ERROR_VARIABLE _err)
  file(WRITE "${BINARY_DIR}/launcher-${_case}.log" "${_out}${_err}")
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "Launcher ${_case} failed: ${_out}${_err}")
  endif()
  file(READ "${_trace}" _actual)
  if(NOT _actual STREQUAL _expected)
    message(FATAL_ERROR "Launcher ${_case}: expected '${_expected}', got '${_actual}'")
  endif()
endforeach()
message(STATUS "PCH on/off, mixed C/C++, consumer and install checks passed")

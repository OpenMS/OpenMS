# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.24)

# cmake -DBINARY_DIR=<scratch> [-DGENERATOR="Unix Makefiles"] -P run.cmake
if(NOT BINARY_DIR)
  message(FATAL_ERROR "Set BINARY_DIR to a scratch directory")
endif()
get_filename_component(BINARY_DIR "${BINARY_DIR}" ABSOLUTE)
file(MAKE_DIRECTORY "${BINARY_DIR}")
if(NOT GENERATOR)
  set(GENERATOR Ninja)
endif()
set(_source "${CMAKE_CURRENT_LIST_DIR}")
set(_configure -G "${GENERATOR}" -DCMAKE_BUILD_TYPE=Release)
foreach(_variable IN ITEMS CMAKE_CXX_COMPILER CMAKE_MAKE_PROGRAM CMAKE_OSX_ARCHITECTURES
    CMAKE_OSX_SYSROOT CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_MSVC_RUNTIME_LIBRARY
    CMAKE_PREFIX_PATH CMAKE_TOOLCHAIN_FILE VCPKG_INSTALLED_DIR VCPKG_TARGET_TRIPLET VCPKG_HOST_TRIPLET)
  if(DEFINED ${_variable} AND NOT "${${_variable}}" STREQUAL "")
    string(REPLACE ";" "\\;" _value "${${_variable}}")
    list(APPEND _configure "-D${_variable}=${_value}")
  endif()
endforeach()
if(GENERATOR_PLATFORM)
  list(APPEND _configure -A "${GENERATOR_PLATFORM}")
endif()
if(GENERATOR_TOOLSET)
  list(APPEND _configure -T "${GENERATOR_TOOLSET}")
endif()

function(run step)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE _result OUTPUT_VARIABLE _out ERROR_VARIABLE _err)
  file(WRITE "${BINARY_DIR}/${step}.log" "${_out}${_err}")
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "${step} failed (${_result}):\n${_out}${_err}")
  endif()
endfunction()

set(_build "${BINARY_DIR}/producer")
set(_prefix "${BINARY_DIR}/prefix")
run(configure "${CMAKE_COMMAND}" -S "${_source}" -B "${_build}" ${_configure}
  -DOPENMS_VERIFY_INTERFACE_HEADER_SETS=OFF)
run(build "${CMAKE_COMMAND}" --build "${_build}" --config Release --parallel 2)
run(install "${CMAKE_COMMAND}" --install "${_build}" --config Release --prefix "${_prefix}")

file(GLOB_RECURSE _excluded "${_prefix}/*private.h" "${_prefix}/*OpenMSTestFramework*"
  "${_prefix}/include/OpenMS/CONCEPT/*")
if(_excluded)
  message(FATAL_ERROR "A default install leaked private or test-framework files")
endif()

# Including API.h also checks all three generated headers through the installed export.
run(consumer-configure "${CMAKE_COMMAND}" -S "${_source}/consumer"
  -B "${BINARY_DIR}/consumer" ${_configure}
  "-DTARGETS_FILE=${_prefix}/lib/cmake/OpenMS/OpenMSTargets.cmake")
run(consumer-build "${CMAKE_COMMAND}" --build "${BINARY_DIR}/consumer" --config Release --parallel 2)

set(_headers_prefix "${BINARY_DIR}/headers-only")
run(install-headers "${CMAKE_COMMAND}" --install "${_build}" --config Release
  --prefix "${_headers_prefix}" --component OpenMSHeaderProbe_headers)
run(install-framework-headers "${CMAKE_COMMAND}" --install "${_build}" --config Release
  --prefix "${_headers_prefix}" --component OpenMSTestFramework_headers)
set(_expected
  include/OpenMS/API.h include/OpenMS/OpenMSHeaderProbeConfig.h
  include/OpenMS/build_config.h include/OpenMS/configured.h
  include/OpenMS/CONCEPT/ClassTest.h include/OpenMS/CONCEPT/ClassTestUtils.h
  include/OpenMS/CONCEPT/FuzzyStringComparator.h include/OpenMS/CONCEPT/MacrosTest.h)
list(SORT _expected)
file(GLOB_RECURSE _installed RELATIVE "${_headers_prefix}" "${_headers_prefix}/*")
if(NOT _installed STREQUAL _expected)
  message(FATAL_ERROR "Unexpected header-component installation: ${_installed}")
endif()

# With verification on, every public header is compiled as its own translation
# unit and the private file set must stay out of it: source/private.h and
# include/OpenMS/private.h are #error headers, so a private header leaking into
# an interface file set fails this build instead of passing unnoticed.
run(verify-configure "${CMAKE_COMMAND}" -S "${_source}" -B "${BINARY_DIR}/verify" ${_configure}
  -DOPENMS_VERIFY_INTERFACE_HEADER_SETS=ON)
run(verify "${CMAKE_COMMAND}" --build "${BINARY_DIR}/verify" --config Release
  --target all_verify_interface_header_sets --parallel 2)

# The configure-time guard must remain active with verification disabled.
execute_process(COMMAND "${CMAKE_COMMAND}"
  -S "${_source}" -B "${BINARY_DIR}/json" ${_configure}
  -DBAD_HEADER=json -DOPENMS_VERIFY_INTERFACE_HEADER_SETS=OFF
  RESULT_VARIABLE _result OUTPUT_VARIABLE _out ERROR_VARIABLE _err)
file(WRITE "${BINARY_DIR}/reject-json.log" "${_out}${_err}")
# CMake wraps diagnostics according to path length; whitespace may include a newline.
if(_result EQUAL 0 OR NOT "${_out}${_err}" MATCHES "includes[ \t\r\n]+nlohmann/json")
  message(FATAL_ERROR "JSON dependency was not rejected as expected:\n${_out}${_err}")
endif()
message(STATUS "Header installation, consumer compatibility, header verification and JSON guard passed")

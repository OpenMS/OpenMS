# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.24)

# cmake -DBINARY_DIR=<scratch> [-DTEST_QT=ON]
#       [-DCONSUMER_CMAKE=<cmake-3.22>] -P run.cmake
if(NOT BINARY_DIR)
  message(FATAL_ERROR "Set BINARY_DIR to a scratch directory")
endif()
get_filename_component(BINARY_DIR "${BINARY_DIR}" ABSOLUTE)
file(MAKE_DIRECTORY "${BINARY_DIR}")
if(NOT GENERATOR)
  set(GENERATOR Ninja)
endif()
if(NOT CONSUMER_CMAKE)
  set(CONSUMER_CMAKE "${CMAKE_COMMAND}")
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

function(expect_failure step pattern)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE _result OUTPUT_VARIABLE _out ERROR_VARIABLE _err)
  file(WRITE "${BINARY_DIR}/${step}.log" "${_out}${_err}")
  if(_result EQUAL 0 OR NOT "${_out}${_err}" MATCHES "${pattern}")
    message(FATAL_ERROR "${step} did not fail as expected:\n${_out}${_err}")
  endif()
endfunction()

function(assert_header prefix relative original)
  if(NOT EXISTS "${prefix}/include/${relative}")
    message(FATAL_ERROR "Missing installed header: ${relative}")
  endif()
  file(SHA256 "${prefix}/include/${relative}" _installed)
  file(SHA256 "${original}" _original)
  if(NOT _installed STREQUAL _original)
    message(FATAL_ERROR "Installed header changed: ${relative}")
  endif()
endfunction()

set(_build "${BINARY_DIR}/producer")
set(_prefix "${BINARY_DIR}/prefix")
run(configure "${CMAKE_COMMAND}" -S "${_source}" -B "${_build}" ${_configure}
  -DOPENMS_VERIFY_INTERFACE_HEADER_SETS=ON "-DTEST_QT=${TEST_QT}")
run(build "${CMAKE_COMMAND}" --build "${_build}" --config Release --parallel 2)
run(verify "${CMAKE_COMMAND}" --build "${_build}" --config Release
  --target all_verify_interface_header_sets --parallel 2)
run(install "${CMAKE_COMMAND}" --install "${_build}" --config Release --prefix "${_prefix}")

assert_header("${_prefix}" OpenMS/API.h "${_source}/include/OpenMS/API.h")
foreach(_header IN ITEMS configured.h build_config.h OpenMSHeaderProbeConfig.h)
  assert_header("${_prefix}" "OpenMS/${_header}" "${_build}/include/OpenMS/${_header}")
endforeach()
if(TEST_QT)
  assert_header("${_prefix}" OpenMS/ProbeObject.h "${_source}/include/OpenMS/ProbeObject.h")
  file(GLOB_RECURSE _moc "${_build}/OpenMSHeaderProbe_autogen/moc_ProbeObject.cpp")
  if(NOT _moc)
    message(FATAL_ERROR "AUTOMOC did not process the public file-set header")
  endif()
endif()
file(GLOB_RECURSE _private "${_prefix}/*private.h")
file(GLOB_RECURSE _framework_library "${_prefix}/*OpenMSTestFramework*")
if(_private OR _framework_library OR EXISTS "${_prefix}/include/OpenMS/CONCEPT/ClassTest.h")
  message(FATAL_ERROR "A default install leaked private or test-framework files")
endif()

# Build-tree and installed exports both support older consumers. Building these
# consumers checks actual include lookup and linkage, not just target properties.
foreach(_kind IN ITEMS build installed)
  if(_kind STREQUAL build)
    set(_targets "${_build}/OpenMSTargets.cmake")
  else()
    set(_targets "${_prefix}/lib/cmake/OpenMS/OpenMSTargets.cmake")
  endif()
  run(consumer-${_kind}-configure "${CONSUMER_CMAKE}" -S "${_source}/consumer"
    -B "${BINARY_DIR}/consumer-${_kind}" ${_configure}
    "-DTARGETS_FILE=${_targets}" "-DTEST_QT=${TEST_QT}")
  run(consumer-${_kind}-build "${CONSUMER_CMAKE}" --build "${BINARY_DIR}/consumer-${_kind}"
    --config Release --parallel 2)
endforeach()

set(_headers_prefix "${BINARY_DIR}/headers-only")
run(install-headers "${CMAKE_COMMAND}" --install "${_build}" --config Release
  --prefix "${_headers_prefix}" --component OpenMSHeaderProbe_headers)
assert_header("${_headers_prefix}" OpenMS/API.h "${_source}/include/OpenMS/API.h")
file(GLOB_RECURSE _libraries "${_headers_prefix}/lib/*")
if(_libraries)
  message(FATAL_ERROR "Installing the header component also installed libraries")
endif()
run(install-framework-headers "${CMAKE_COMMAND}" --install "${_build}" --config Release
  --prefix "${_headers_prefix}" --component OpenMSTestFramework_headers)
get_filename_component(_repo "${_source}/../../.." ABSOLUTE)
foreach(_header IN ITEMS ClassTest.h ClassTestUtils.h FuzzyStringComparator.h MacrosTest.h)
  assert_header("${_headers_prefix}" "OpenMS/CONCEPT/${_header}"
    "${_repo}/src/testframework/include/OpenMS/CONCEPT/${_header}")
endforeach()

# The configure-time guard must remain active with verification disabled.
expect_failure(reject-json "includes nlohmann/json" "${CMAKE_COMMAND}"
  -S "${_source}" -B "${BINARY_DIR}/json" ${_configure}
  -DBAD_HEADER=json -DOPENMS_VERIFY_INTERFACE_HEADER_SETS=OFF)

# A normal build does not compile the verification units; an explicit check
# must fail if a public header needs a PRIVATE dependency's include directory.
set(_leak_build "${BINARY_DIR}/private-dependency")
run(leak-configure "${CMAKE_COMMAND}" -S "${_source}" -B "${_leak_build}" ${_configure}
  -DBAD_HEADER=private_dependency -DOPENMS_VERIFY_INTERFACE_HEADER_SETS=ON)
run(leak-default-build "${CMAKE_COMMAND}" --build "${_leak_build}" --config Release --parallel 2)
expect_failure(reject-private-dependency "private_dependency.h" "${CMAKE_COMMAND}"
  --build "${_leak_build}" --config Release --target all_verify_interface_header_sets --parallel 2)
message(STATUS "Header file-set installation, consumers, verification and exclusions passed")

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.24)
project(OpenMSHeaderSetsTest CXX)

# Exercise the production library/install helpers without building OpenMS or its
# dependencies. The test-framework target below is the real, standalone library.
get_filename_component(OPENMS_HOST_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(OPENMS_HOST_BINARY_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
list(APPEND CMAKE_MODULE_PATH "${OPENMS_HOST_DIRECTORY}/cmake/Modules")
set(INSTALL_LIB_DIR lib)
set(INSTALL_INCLUDE_DIR include)
set(INSTALL_CMAKE_DIR lib/cmake/OpenMS)
set(BUILD_SHARED_LIBS ON)
include("${OPENMS_HOST_DIRECTORY}/cmake/install_macros.cmake")
include("${OPENMS_HOST_DIRECTORY}/cmake/export_macros.cmake")
include("${OPENMS_HOST_DIRECTORY}/cmake/add_library_macros.cmake")
function(openms_add_library_compiler_flags target)
endfunction()
function(openms_doc_path path)
endfunction()

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/include/OpenMS/configured.h" "#pragma once\n#define PROBE_VALUE 42\n")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/build_config.h.in" "#pragma once\n")
set(_build_header "${CMAKE_CURRENT_BINARY_DIR}/include/OpenMS/build_config.h")
add_custom_command(OUTPUT "${_build_header}"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CMAKE_CURRENT_BINARY_DIR}/build_config.h.in" "${_build_header}"
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/build_config.h.in" VERBATIM)
set(_headers include/OpenMS/API.h
  "${CMAKE_CURRENT_BINARY_DIR}/include/OpenMS/configured.h" "${_build_header}")
if(BAD_HEADER STREQUAL "json")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/include/OpenMS/bad.h" "#include <nlohmann/json.hpp>\n")
  list(APPEND _headers "${CMAKE_CURRENT_BINARY_DIR}/include/OpenMS/bad.h")
elseif(BAD_HEADER STREQUAL "private_dependency")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/include/OpenMS/bad.h" "#include <private_dependency.h>\n")
  list(APPEND _headers "${CMAKE_CURRENT_BINARY_DIR}/include/OpenMS/bad.h")
endif()

openms_add_library(TARGET_NAME OpenMSHeaderProbe
  SOURCE_FILES probe.cpp
  HEADER_FILES ${_headers}
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}/include"
  PRIVATE_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/private_dependency"
  DLL_EXPORT_PATH "OpenMS/")
target_sources(OpenMSHeaderProbe PRIVATE FILE_SET private_headers TYPE HEADERS
  BASE_DIRS source include
  FILES source/private.h include/OpenMS/private.h)
if(NOT OpenMSHeaderProbe_INCLUDE_DIRECTORIES STREQUAL
    "${CMAKE_CURRENT_SOURCE_DIR}/include;${CMAKE_CURRENT_BINARY_DIR}/include")
  message(FATAL_ERROR "The class-test include-directory cache lost its roots")
endif()

if(TEST_QT)
  find_package(Qt6 REQUIRED COMPONENTS Core)
  target_sources(OpenMSHeaderProbe PUBLIC FILE_SET HEADERS FILES include/OpenMS/ProbeObject.h)
  target_link_libraries(OpenMSHeaderProbe PUBLIC Qt6::Core)
  set_target_properties(OpenMSHeaderProbe PROPERTIES AUTOMOC ON)
  # A missing moc output leaves this reference unresolved when the consumer links.
  target_sources(OpenMSHeaderProbe PRIVATE probe_object.cpp)
endif()

add_subdirectory("${OPENMS_HOST_DIRECTORY}/src/testframework" testframework)
install_export_targets()
export(EXPORT OpenMSTargets NAMESPACE OpenMS:: FILE "${CMAKE_CURRENT_BINARY_DIR}/OpenMSTargets.cmake")

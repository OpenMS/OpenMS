# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# $Maintainer: Karnajeet Gosavi $

# Called after project() has loaded the compiler's module-discovery result.
function(openms_validate_import_std)
  if(NOT OPENMS_IMPORT_STD)
    return()
  endif()
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(FATAL_ERROR "OPENMS_IMPORT_STD is an MSVC-only experiment; use cl.exe or disable the option.")
  endif()
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.5)
    # CMake 4.5 replaced the supported-standard list with module metadata.
    if(EXISTS "${CMAKE_CXX_STDLIB_MODULES_JSON}")
      return()
    endif()
    set(_reason "${CMAKE_CXX_COMPILER_IMPORT_STD_ERROR_MESSAGE}")
  else()
    if("23" IN_LIST CMAKE_CXX_COMPILER_IMPORT_STD)
      return()
    endif()
    set(_reason "${CMAKE_CXX23_COMPILER_IMPORT_STD_NOT_FOUND_MESSAGE}")
  endif()
  message(FATAL_ERROR
    "OPENMS_IMPORT_STD: CMake did not detect C++23 standard-library module support.\n"
    "Compiler detection: ${_reason}\n"
    "When enabling OPENMS_IMPORT_STD in an existing build, configure a fresh build directory "
    "to refresh compiler/module detection.")
endfunction()

option(OPENMS_IMPORT_STD "Experiment with MSVC C++23 import std in selected OpenMS sources." OFF)
if(NOT OPENMS_IMPORT_STD)
  return()
endif()

if(CMAKE_VERSION VERSION_LESS 3.30)
  message(FATAL_ERROR "OPENMS_IMPORT_STD requires CMake 3.30 or newer.")
endif()
if(NOT CMAKE_GENERATOR MATCHES "^Ninja")
  message(FATAL_ERROR "OPENMS_IMPORT_STD requires Ninja or Ninja Multi-Config; Visual Studio generators cannot build imported standard library modules.")
endif()

# Values from Help/dev/experimental.rst in the corresponding CMake releases.
# This must run before project() detects the C++ toolchain. Some UUID changes
# were backported, so the patch versions matter as well as the minor versions.
if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.5)
  # CMake 4.5 no longer requires the experimental gate.
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 4.4)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "f35a9ac6-8463-4d38-8eec-5d6008153e7d")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 4.3)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "451f2fe2-a8a2-47c3-bc32-94786d8fc91b")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.3 OR
       (CMAKE_VERSION VERSION_GREATER_EQUAL 3.31.8 AND CMAKE_VERSION VERSION_LESS 4.0))
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "a9e1cf81-9932-4810-974b-6eccaf14e457")
else()
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")
endif()

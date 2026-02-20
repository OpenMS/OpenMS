# cmake/import_std_setup.cmake
#
# Maps CMake versions to the experimental UUID needed for import std support.
# Must be included BEFORE project().
#
# See: https://www.kitware.com/import-std-in-cmake-3-30/
# See: https://cmake.org/cmake/help/latest/prop_tgt/CXX_MODULE_STD.html

# Known UUIDs per CMake version range.
# Check Help/dev/experimental.rst in CMake source when adding new versions.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "4.1")
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0" AND CMAKE_VERSION VERSION_LESS "4.1")
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "a9e1cf81-9932-4810-974b-6eccaf14e457")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30" AND CMAKE_VERSION VERSION_LESS "4.0")
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")
else()
  message(FATAL_ERROR
    "No known CMAKE_EXPERIMENTAL_CXX_IMPORT_STD UUID for CMake ${CMAKE_VERSION}. "
    "Check Help/dev/experimental.rst in the CMake source for your version and update "
    "cmake/import_std_setup.cmake accordingly.")
endif()

message(STATUS "import std: enabled (CMake ${CMAKE_VERSION}, UUID ${CMAKE_EXPERIMENTAL_CXX_IMPORT_STD})")

# Enable import std for all targets
set(CMAKE_CXX_MODULE_STD ON)

# Enable C++ module source scanning for all targets.
# Required because subdirectories use cmake_minimum_required(VERSION 3.15)
# which sets CMP0155 to OLD (disabling scanning). This global setting
# ensures all targets get module dependency scanning regardless of policy.
set(CMAKE_CXX_SCAN_FOR_MODULES ON)

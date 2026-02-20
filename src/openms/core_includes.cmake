# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# core_includes.cmake — aggregates sources for OpenMS_Core library
# Core modules: INTERFACES, CONCEPT, SYSTEM, MATH, DATASTRUCTURES, METADATA, KERNEL, CHEMISTRY, IONMOBILITY
#
# Note: ${CMAKE_CURRENT_LIST_DIR} is needed because this file is included from
# a subdirectory CMakeLists.txt whose CMAKE_CURRENT_SOURCE_DIR differs.

set(OpenMS_sources CACHE INTERNAL "")

## Source files — order respects dependency hierarchy
include(${CMAKE_CURRENT_LIST_DIR}/source/INTERFACES_IMPL/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/CONCEPT/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/SYSTEM/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/MATH/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/MATH/STATISTICS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/MATH/MISC/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/DATASTRUCTURES/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/METADATA/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/METADATA/ID/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/KERNEL/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/CHEMISTRY/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/CHEMISTRY/ISOTOPEDISTRIBUTION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/CHEMISTRY/MASSDECOMPOSITION/IMS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/CHEMISTRY/MASSDECOMPOSITION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/IONMOBILITY/sources.cmake)

set(OpenMS_Core_sources ${OpenMS_sources} CACHE INTERNAL "OpenMS Core source files")
set(OpenMS_sources CACHE INTERNAL "")

## Header files
set(OpenMS_sources_h CACHE INTERNAL "")

include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/INTERFACES/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/CONCEPT/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/SYSTEM/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/MATH/MISC/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/MATH/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/MATH/STATISTICS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/DATASTRUCTURES/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/DATASTRUCTURES/Utils/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/METADATA/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/METADATA/ID/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/KERNEL/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/CHEMISTRY/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/IONMOBILITY/sources.cmake)

set(OpenMS_Core_sources_h ${OpenMS_sources_h} CACHE INTERNAL "OpenMS Core header files")
set(OpenMS_sources_h CACHE INTERNAL "")

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

# Remove source files that have IO/Algo back-dependencies and belong in those layers.
# These files are leaf nodes in Core's dependency graph (no other Core .cpp depends on them)
# but reference FORMAT types. They are compiled into OpenMS_IO or OpenMS_Algo instead.
set(_core_to_io_sources
  source/CONCEPT/ClassTest.cpp        # test framework validation uses FORMAT file classes
  source/CONCEPT/FuzzyStringComparator.cpp  # uses TextFile
  source/KERNEL/OnDiscMSExperiment.cpp      # uses FileHandler, PeakFileOptions, IndexedMzMLHandler
  source/METADATA/ID/IdentificationDataConverter.cpp  # uses MzTab types
  source/METADATA/ExperimentalDesign.cpp    # uses FileHandler, TextFile
)
set(_core_to_algo_sources
  source/MATH/STATISTICS/PosteriorErrorProbabilityModel.cpp  # uses TextFile for gnuplot output
)
list(REMOVE_ITEM OpenMS_sources ${_core_to_io_sources} ${_core_to_algo_sources})

# Capture accumulated sources into Core variable and reset accumulator
set(OpenMS_Core_sources ${OpenMS_sources} CACHE INTERNAL "OpenMS Core source files")
set(OpenMS_Core_moved_to_io ${_core_to_io_sources} CACHE INTERNAL "Core sources moved to IO")
set(OpenMS_Core_moved_to_algo ${_core_to_algo_sources} CACHE INTERNAL "Core sources moved to Algo")
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

# Remove headers for files moved to IO/Algo (parallel to source file moves above)
set(_core_to_io_headers
  include/OpenMS/CONCEPT/ClassTest.h
  include/OpenMS/CONCEPT/FuzzyStringComparator.h
  include/OpenMS/KERNEL/OnDiscMSExperiment.h
  include/OpenMS/METADATA/ID/IdentificationDataConverter.h
  include/OpenMS/METADATA/ExperimentalDesign.h
)
set(_core_to_algo_headers
  include/OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h
)
list(REMOVE_ITEM OpenMS_sources_h ${_core_to_io_headers} ${_core_to_algo_headers})

set(OpenMS_Core_sources_h ${OpenMS_sources_h} CACHE INTERNAL "OpenMS Core header files")
set(OpenMS_Core_moved_to_io_h ${_core_to_io_headers} CACHE INTERNAL "Core headers moved to IO")
set(OpenMS_Core_moved_to_algo_h ${_core_to_algo_headers} CACHE INTERNAL "Core headers moved to Algo")
set(OpenMS_sources_h CACHE INTERNAL "")

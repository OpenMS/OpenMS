# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# io_includes.cmake — aggregates sources for OpenMS_IO library
# IO modules: FORMAT and all subdirectories
#
# Note: ${CMAKE_CURRENT_LIST_DIR} is needed because this file is included from
# a subdirectory CMakeLists.txt whose CMAKE_CURRENT_SOURCE_DIR differs.

set(OpenMS_sources CACHE INTERNAL "")

## Source files
include(${CMAKE_CURRENT_LIST_DIR}/source/FORMAT/DATAACCESS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/FORMAT/HANDLERS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/FORMAT/MSNUMPRESS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/FORMAT/VALIDATORS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/FORMAT/OPTIONS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/FORMAT/sources.cmake)

# Add Core source files that were moved here due to IO back-dependencies
list(APPEND OpenMS_sources ${OpenMS_Core_moved_to_io})

# Remove source files that have Algo back-dependencies and belong in the Algo layer.
# These FORMAT files depend on ANALYSIS/TOPDOWN types (PeakGroup, DeconvolvedSpectrum, etc.)
set(_io_to_algo_sources
  source/FORMAT/FLASHDeconvFeatureFile.cpp   # uses PeakGroup, DeconvolvedSpectrum, FLASHHelperClasses
  source/FORMAT/FLASHDeconvSpectrumFile.cpp  # uses PeakGroup, DeconvolvedSpectrum, FLASHHelperClasses
)
list(REMOVE_ITEM OpenMS_sources ${_io_to_algo_sources})

set(OpenMS_IO_sources ${OpenMS_sources} CACHE INTERNAL "OpenMS IO source files")
set(OpenMS_IO_moved_to_algo ${_io_to_algo_sources} CACHE INTERNAL "IO sources moved to Algo")
set(OpenMS_sources CACHE INTERNAL "")

## Header files
set(OpenMS_sources_h CACHE INTERNAL "")

include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/DATAACCESS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/HANDLERS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/MSNUMPRESS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/VALIDATORS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/OPTIONS/sources.cmake)

# Add Core headers that were moved here (parallel to source moves)
list(APPEND OpenMS_sources_h ${OpenMS_Core_moved_to_io_h})

# Remove headers for files compiled in Algo (parallel to source file moves above)
set(_io_to_algo_headers
  include/OpenMS/FORMAT/FLASHDeconvFeatureFile.h
  include/OpenMS/FORMAT/FLASHDeconvSpectrumFile.h
)
list(REMOVE_ITEM OpenMS_sources_h ${_io_to_algo_headers})

set(OpenMS_IO_sources_h ${OpenMS_sources_h} CACHE INTERNAL "OpenMS IO header files")
set(OpenMS_IO_moved_to_algo_h ${_io_to_algo_headers} CACHE INTERNAL "IO headers moved to Algo")
set(OpenMS_sources_h CACHE INTERNAL "")

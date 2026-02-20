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

# Core source files compiled here because they depend on IO (FORMAT) types
# or IO-tier libraries. The files remain in their original directories but are
# not listed in Core's sources.cmake — they appear only here.
list(APPEND OpenMS_sources
  source/CONCEPT/ClassTest.cpp                         # uses FORMAT file classes
  source/CONCEPT/FuzzyStringComparator.cpp             # uses TextFile
  source/CONCEPT/Init.cpp                              # Xerces-C init (only IO uses Xerces)
  source/METADATA/ID/IdentificationDataConverter.cpp   # uses MzTab types
  source/METADATA/ExperimentalDesign.cpp               # uses FileHandler, TextFile
  source/METADATA/SpectrumMetaDataLookup.cpp           # uses FileHandler
)

set(OpenMS_IO_sources ${OpenMS_sources} CACHE INTERNAL "OpenMS IO source files")
set(OpenMS_sources CACHE INTERNAL "")

## Header files
set(OpenMS_sources_h CACHE INTERNAL "")

include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/DATAACCESS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/HANDLERS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/MSNUMPRESS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/VALIDATORS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FORMAT/OPTIONS/sources.cmake)

# Core headers for the source files above (parallel to source moves)
list(APPEND OpenMS_sources_h
  include/OpenMS/CONCEPT/ClassTest.h
  include/OpenMS/CONCEPT/FuzzyStringComparator.h
  include/OpenMS/CONCEPT/Init.h
  include/OpenMS/METADATA/ID/IdentificationDataConverter.h
  include/OpenMS/METADATA/ExperimentalDesign.h
  include/OpenMS/METADATA/SpectrumMetaDataLookup.h
)

set(OpenMS_IO_sources_h ${OpenMS_sources_h} CACHE INTERNAL "OpenMS IO header files")
set(OpenMS_sources_h CACHE INTERNAL "")

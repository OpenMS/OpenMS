# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


project("OpenMS-IsoSpec")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)


set(_ISOSPEC_SOURCES
  cwrapper.cpp
  allocator.cpp  
  dirtyAllocator.cpp 
  isoSpec++.cpp 
  isoMath.cpp 
  marginalTrek++.cpp 
  operators.cpp
  element_tables.cpp
  misc.cpp
  fasta.cpp
)

set(ISOSPEC_SOURCES )
foreach(_file ${_ISOSPEC_SOURCES})
  # append source files
  set(
    ISOSPEC_SOURCES
    ${ISOSPEC_SOURCES}
    ${PROJECT_SOURCE_DIR}/IsoSpec/${_file}
  )
endforeach()

# Create IsoSpec as a proper library target
add_library(IsoSpec STATIC ${ISOSPEC_SOURCES})

# Set include directories for the library
target_include_directories(IsoSpec SYSTEM PUBLIC 
                          "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/IsoSpec>"
                          "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty/IsoSpec>")

# Apply compiler flags
if (NOT WIN32)
  target_compile_options(IsoSpec PRIVATE -fPIC)
endif()

# Install the library
install_library(IsoSpec)
openms_register_export_target(IsoSpec)

# Install headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/IsoSpec/
  DESTINATION ${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty/IsoSpec
  COMPONENT thirdparty_headers
  FILES_MATCHING PATTERN "*.h"
  PATTERN "*.cpp" EXCLUDE
)

get_filename_component(ISOSPEC_FOLDER "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(APPEND ${TPLICENSEFILE} "\nIsoSpec:\n")
cat(${ISOSPEC_FOLDER}/LICENSE ${TPLICENSEFILE})

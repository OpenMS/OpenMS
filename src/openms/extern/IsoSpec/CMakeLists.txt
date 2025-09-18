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
)

# export include directory
set(ISOSPEC_INCLUDE_DIRECTORY
  ${PROJECT_SOURCE_DIR}
  CACHE
  INTERNAL "isospec include directories" FORCE
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

# TODO Instead of including including the complete sources, make a library,
#  link, install and export
# Commented out until used
#add_library(IsoSpec ${ISOSPEC_SOURCES})
#target_include_directories(IsoSpec SYSTEM INTERFACE ${ISOSPEC_INCLUDE_DIRECTORY})

get_filename_component(ISOSPEC_FOLDER "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(APPEND ${TPLICENSEFILE} "\nIsoSpec:\n")
cat(${ISOSPEC_FOLDER}/LICENSE ${TPLICENSEFILE})

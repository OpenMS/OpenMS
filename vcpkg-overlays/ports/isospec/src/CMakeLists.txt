# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


project("IsoSpec")
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

# Add mman.cpp on Windows for memory mapping support
if(WIN32)
  list(APPEND _ISOSPEC_SOURCES mman.cpp)
endif()

set(ISOSPEC_SOURCES )
foreach(_file ${_ISOSPEC_SOURCES})
  # append source files
  set(
    ISOSPEC_SOURCES
    ${ISOSPEC_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/IsoSpec/${_file}
  )
endforeach()

# Create IsoSpec as a proper library target
add_library(IsoSpec STATIC ${ISOSPEC_SOURCES})  

# Require C++11 for IsoSpec
target_compile_features(IsoSpec PUBLIC cxx_std_11)

target_compile_definitions(IsoSpec PUBLIC
  ISOSPEC_GOT_SYSTEM_MMAN=false
  ISOSPEC_GOT_MMAN=false
  ISOSPEC_BUILDING_OPENMS=true
)

# Set include directories for the library
target_include_directories(IsoSpec SYSTEM PUBLIC 
                          "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
                          "$<INSTALL_INTERFACE:include>")

# Generate export header for Windows DLL symbols
include(GenerateExportHeader)
generate_export_header(IsoSpec
                      BASE_NAME OPENMS
                      EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/IsoSpec/IsoSpecExport.h)

# Add the binary directory to include path so the generated header can be found
target_include_directories(IsoSpec PUBLIC
                          "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>")

# Install the generated export header
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/IsoSpec/IsoSpecExport.h
        DESTINATION include/IsoSpec)

# Objects are linked into a shared library, so they must be position independent.
set_target_properties(IsoSpec PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS IsoSpec
  EXPORT IsoSpecTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(EXPORT IsoSpecTargets
  FILE IsoSpecConfig.cmake
  DESTINATION share/isospec
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/IsoSpec/
  DESTINATION include/IsoSpec
  FILES_MATCHING PATTERN "*.h"
  PATTERN "*.cpp" EXCLUDE
)
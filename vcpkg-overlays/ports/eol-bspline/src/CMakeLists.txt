# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


project("eol-bspline")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(_EOL_BSPLINE_SOURCES
BSplineLib.cpp
)

set(EOL_BSPLINE_SOURCES )
foreach(_file ${_EOL_BSPLINE_SOURCES})
  # append source files
  set(
    EOL_BSPLINE_SOURCES
    ${EOL_BSPLINE_SOURCES}
    ${PROJECT_SOURCE_DIR}/BSpline/${_file}
  )
endforeach()

# Create eol-bspline as a proper library target
add_library(eol-bspline STATIC ${EOL_BSPLINE_SOURCES})

# Set include directories for the library
target_include_directories(eol-bspline SYSTEM PUBLIC 
                          "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
                          "$<INSTALL_INTERFACE:include>")

# Apply compiler flags
if (NOT WIN32)
  target_compile_options(eol-bspline PRIVATE -fPIC)
endif()

install(TARGETS eol-bspline
  EXPORT eol-bsplineTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(EXPORT eol-bsplineTargets
  FILE eol-bsplineConfig.cmake
  DESTINATION share/eol-bspline
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/BSpline/
  DESTINATION include/BSpline
  FILES_MATCHING PATTERN "*.h"
  PATTERN "*.cpp" EXCLUDE
)

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


project("OpenMS-eol-bsplines")
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
                          "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty>")

# Apply compiler flags
if (NOT WIN32)
  target_compile_options(eol-bspline PRIVATE -fPIC)
endif()

# Install the library
install_library(eol-bspline)
openms_register_export_target(eol-bspline)

# Install headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/BSpline/
  DESTINATION ${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty/BSpline
  COMPONENT thirdparty_headers
  FILES_MATCHING PATTERN "*.h"
  PATTERN "*.cpp" EXCLUDE
)

file(APPEND ${TPLICENSEFILE} "\neol-bspline:\n")
cat(${PROJECT_SOURCE_DIR}/LICENSE ${TPLICENSEFILE})

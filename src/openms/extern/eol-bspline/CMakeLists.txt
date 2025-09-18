# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


project("OpenMS-eol-bsplines")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# export include directory
set(EOL_BSPLINE_INCLUDE_DIRECTORY
  ${PROJECT_SOURCE_DIR}
  CACHE
  INTERNAL "eol-bspline include directories" FORCE
)

set(_EOL_BSPLINE_SOURCES
BSpline.cpp
BSpline.h
BSplineBase.cpp
BSplineBase.h
BSplineLib.cpp
BandedMatrix.h
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

file(APPEND ${TPLICENSEFILE} "\neol-bspline:\n")
cat(${PROJECT_SOURCE_DIR}/LICENSE ${TPLICENSEFILE})

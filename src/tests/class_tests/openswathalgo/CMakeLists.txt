# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Witold Wolski $
# $Authors: Witold Wolski, Stephan Aiche $
# --------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_class_tests_openswathalgo")

# --------------------------------------------------------------------------
# add include directories
include_directories(SYSTEM ${OpenMS_INCLUDE_DIRECTORIES})

# --------------------------------------------------------------------------
# list of available tests
set(openswath_algo_tests
  Scoring_test
  Datastructures_test
  TestConvert
  DiaHelpers_test
  SwathMap_test
)

#------------------------------------------------------------------------------
# QT dependencies
#find_package(Qt5 COMPONENTS Core Network REQUIRED)
#if (NOT Qt5Network_FOUND)
#  message(STATUS "QtNetwork module not found!")
#  message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
#endif()

# --------------------------------------------------------------------------
# add targets for the executables
foreach(i ${openswath_algo_tests})
  add_executable(${i} ${i}.cpp)
  # OpenMS is listed explicitly and *after* the framework: the framework still uses
  # libOpenMS' low-level utilities (StringUtils, Exception, ...), and these test
  # sources themselves reference no OpenMS symbols, so with the linker's default
  # --as-needed libOpenMS.so has to follow the archive whose references it satisfies.
  # Once libOpenMS is split up, this becomes a dependency on the future core library.
  target_link_libraries(${i} OpenSwathAlgo OpenMSTestFramework OpenMS)
  add_test(${i} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i})
endforeach(i)

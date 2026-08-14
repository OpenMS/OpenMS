# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_class_tests")

#------------------------------------------------------------------------------
# set new CMAKE_RUNTIME_OUTPUT_DIRECTORY for tests and remember old setting
set(_TMP_CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#------------------------------------------------------------------------------
# libOpenMS-specific test-framework hooks (UniqueIdGenerator seeding, XML
# validation of tmp files), linked into the openms and openms_gui class tests
add_subdirectory(support)

#------------------------------------------------------------------------------
# add OpenSWATHAlgo specific tests
add_subdirectory(openswathalgo)

#------------------------------------------------------------------------------
# add OpenMS specific tests
add_subdirectory(openms)

if(WITH_GUI)
  #------------------------------------------------------------------------------
  # add OpenMS_GUI specific tests
  add_subdirectory(openms_gui)
endif()

#------------------------------------------------------------------------------
# restore old CMAKE_RUNTIME_OUTPUT_DIRECTORY
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${_TMP_CMAKE_RUNTIME_OUTPUT_DIRECTORY})

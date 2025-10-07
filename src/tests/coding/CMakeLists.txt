# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_coding_tests")

# --------------------------------------------------------------------------
# OpenMS internal style tests
if (NOT WIN32)
  # Check for unintended copies in C++11 for loops
  add_test(
    NAME potential_for_loop_copies
    COMMAND sh -c "grep -R for --exclude-dir=contrib --exclude-dir=cppcheck --exclude-dir=thirdparty ${OPENMS_HOST_DIRECTORY} \
    | grep auto | grep -v begin | grep -v areaBegin | grep -v MZBegin | grep -v OMS_CODING_TEST_EXCLUDE | grep ':' \
    | grep -v '&' | grep -v '//' | grep 'for[ ]*(' || echo 'Success' "
  )

  set_tests_properties(potential_for_loop_copies PROPERTIES
    PASS_REGULAR_EXPRESSION "Success"
  )
endif()
# --------------------------------------------------------------------------
# cppcheck tests
find_package(cppcheck)
if( CPPCHECK_FOUND )
include( ${PROJECT_SOURCE_DIR}/cppcheck.cmake )
endif(CPPCHECK_FOUND)

# --------------------------------------------------------------------------
# cpplint testing
find_program(PYTHON_EXECUTABLE
         python
         DOC "python executable used to perform coding convention test.")
if("${PYTHON_EXECUTABLE}" STREQUAL "PYTHON_EXECUTABLE-NOTFOUND")
  message(STATUS "Info: The program python could not be found. Coding convention check will not be available! Add the location of python(.exe) to your PATH environment variable.")
else()
  include( ${PROJECT_SOURCE_DIR}/createcpplinttests.cmake )
endif("${PYTHON_EXECUTABLE}" STREQUAL "PYTHON_EXECUTABLE-NOTFOUND")

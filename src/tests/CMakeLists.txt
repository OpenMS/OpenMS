# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche, Chris Bielow $
# --------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_tests")

option(ENABLE_STYLE_TESTING "Enables checking of code convention violations (cpplint) and static code analysis (cppchecker). Note that this will disable the regular test system." OFF)
option(ENABLE_TOPP_TESTING "Enables tests for TOPP. Should be disabled only on time constraints (e.g. chunking during continuous integration)." ON)
option(ENABLE_CLASS_TESTING "Enables tests for library classes. Should be disabled only on time constraints (e.g. chunking during continuous integration)." ON)
option(ENABLE_PIPELINE_TESTING "Enables the additional testing of various TOPPAS pipelines when 'make test' is called." ON)

#------------------------------------------------------------------------------
# we only test if we have no package target
if("${PACKAGE_TYPE}" STREQUAL "none")
  if(ENABLE_STYLE_TESTING)
    add_subdirectory(coding)
  else()
    ## configure the regular class ..
    if(ENABLE_CLASS_TESTING)
      add_subdirectory(class_tests)
    endif()  
    # TOPP Test project ..
    if(ENABLE_TOPP_TESTING AND BUILD_TOPP_TOOLS)
      add_subdirectory(topp)
    endif()
    # and if requested also TOPPAS pipeline tests
    if(ENABLE_PIPELINE_TESTING)
      add_subdirectory(toppas)
    endif()
  endif(ENABLE_STYLE_TESTING)
elseif( ENABLE_STYLE_TESTING OR ENABLE_CLASS_TESTING OR ENABLE_TOPP_TESTING OR ENABLE_PIPELINE_TESTING)
  message(WARNING "PACKAGE_TYPE specified along with enabled testing. Tests are only enabled when PACKAGE_TYPE is none. Ignoring tests.")
endif("${PACKAGE_TYPE}" STREQUAL "none")

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
set(_openms_topp_testing_default ON)
if(DEFINED BUILD_TOPP_TOOLS)
  set(_openms_topp_testing_default ${BUILD_TOPP_TOOLS})
endif()
option(ENABLE_TOPP_TESTING "Enables tests for TOPP. Should be disabled only on time constraints (e.g. chunking during continuous integration)." ${_openms_topp_testing_default})
option(ENABLE_CLASS_TESTING "Enables tests for library classes. Should be disabled only on time constraints (e.g. chunking during continuous integration)." ON)
option(ENABLE_PIPELINE_TESTING "Enables the additional testing of various TOPPAS pipelines when 'make test' is called." ${_openms_topp_testing_default})
unset(_openms_topp_testing_default)

if(DEFINED BUILD_TOPP_TOOLS
   AND NOT BUILD_TOPP_TOOLS
   AND (ENABLE_TOPP_TESTING OR ENABLE_PIPELINE_TESTING))
  message(FATAL_ERROR
    "ENABLE_TOPP_TESTING and ENABLE_PIPELINE_TESTING require BUILD_TOPP_TOOLS to be ON.")
endif()
option(ENABLE_BENCHMARK_TESTING "Enables benchmark tests for large datasets (Issue #8788). Requires downloading large external datasets; meant to run on a schedule (e.g. nightly), not on every PR." OFF)

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
    if(ENABLE_TOPP_TESTING)
      add_subdirectory(topp)
    endif()
    # and if requested also TOPPAS pipeline tests
    if(ENABLE_PIPELINE_TESTING)
      add_subdirectory(toppas)
    endif()
    # benchmark tests for large datasets (Issue #8788) - opt-in, downloads large external data
    if(ENABLE_BENCHMARK_TESTING)
      add_subdirectory(benchmarking)
    endif()
  endif(ENABLE_STYLE_TESTING)
elseif( ENABLE_STYLE_TESTING OR ENABLE_CLASS_TESTING OR ENABLE_TOPP_TESTING OR ENABLE_PIPELINE_TESTING OR ENABLE_BENCHMARK_TESTING)
  message(WARNING "PACKAGE_TYPE specified along with enabled testing. Tests are only enabled when PACKAGE_TYPE is none. Ignoring tests.")
endif("${PACKAGE_TYPE}" STREQUAL "none")

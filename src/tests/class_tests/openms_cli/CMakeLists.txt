# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

project("OpenMS_class_tests_openms_cli")

#------------------------------------------------------------------------------
# Class tests of the TOPP tool framework (libOpenMS_CLI). They use the fixtures of
# the core class tests, so the generated test_config.h points at the same data
# directory, and they register libOpenMS behavior with the test framework through
# the same OpenMSTestSupport target.
set(CF_OPENMS_TEST_DATA_PATH "${PROJECT_SOURCE_DIR}/../openms/data/")
set(CONFIGURED_TEST_CONFIG_H "${PROJECT_BINARY_DIR}/include/OpenMS/test_config.h")
configure_file(${PROJECT_SOURCE_DIR}/../openms/include/OpenMS/test_config.h.in ${CONFIGURED_TEST_CONFIG_H})

set(cli_executables_list
  INIUpdater_test
  #MapAlignerBase_test
  ParameterInformation_test
  SearchEngineBase_test
  TOPPBase_test
  TOPPExternalToolBase_test
  ToolHandler_test
  TransformationModelDefaults_test
)

include_directories(${PROJECT_BINARY_DIR}/include ${PROJECT_SOURCE_DIR}/../openms/include)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

#------------------------------------------------------------------------------
# disable optimization for tests for gcc like compilers
if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
  set(_TMP_CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
  set(CMAKE_CXX_FLAGS_RELEASE "-O0")
endif()

foreach(_class_test ${cli_executables_list})
  add_executable(${_class_test} source/${_class_test}.cpp)
  target_link_libraries(${_class_test} OpenMSTestSupport ${OpenMS_CLI_LIBRARIES})
  openms_add_executable_compiler_flags(${_class_test})
  add_test(${_class_test} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_class_test})
  # only add OPENMP flags to gcc linker (except Mac OS X, due to compiler bug
  # see https://sourceforge.net/apps/trac/open-ms/ticket/280 for details)
  if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${_class_test} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endforeach(_class_test)

#------------------------------------------------------------------------------
# restore old compiler flags
if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
  set(CMAKE_CXX_FLAGS_RELEASE "${_TMP_CMAKE_CXX_FLAGS_RELEASE}")
endif()

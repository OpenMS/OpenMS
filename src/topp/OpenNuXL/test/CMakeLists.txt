# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

# Class tests for the NuXL classes that live with the OpenNuXL tool. The class-test
# project adds this directory (see src/tests/class_tests/openms/CMakeLists.txt), so the
# targets below inherit its include directories -- OpenMS/test_config.h among them -- and
# its -O0 test build flags, and are built and run exactly like any other class test.
#
# Each entry lists the tool-local sources its test has to link: the class it covers plus
# the NuXL classes that class calls into.
set(_nuxl_NuXLFragmentAdductDefinition_sources NuXLFragmentAdductDefinition.cpp)
set(_nuxl_NuXLModificationsGenerator_sources NuXLModificationsGenerator.cpp)
set(_nuxl_NuXLParameterParsing_sources NuXLParameterParsing.cpp NuXLFragmentAdductDefinition.cpp
                                       NuXLModificationsGenerator.cpp NuXLPresets.cpp)
set(_nuxl_NuXLReport_sources NuXLReport.cpp NuXLMarkerIonExtractor.cpp NuXLModificationsGenerator.cpp)

foreach(_nuxl_test NuXLFragmentAdductDefinition NuXLModificationsGenerator NuXLParameterParsing NuXLReport)
  set(_nuxl_test_sources)
  foreach(_nuxl_source ${_nuxl_${_nuxl_test}_sources})
    list(APPEND _nuxl_test_sources ${CMAKE_CURRENT_SOURCE_DIR}/../source/${_nuxl_source})
  endforeach()
  # OpenMSTestSupport.cpp registers libOpenMS behavior with the std-only test framework;
  # every class test needs it (see the file's header comment).
  add_executable(${_nuxl_test}_test
    ${_nuxl_test}_test.cpp
    ${_nuxl_test_sources}
    ${OPENMS_CLASS_TESTS_SUPPORT_SOURCE})
  target_include_directories(${_nuxl_test}_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
  target_link_libraries(${_nuxl_test}_test OpenMSTestFramework ${OpenMS_LIBRARIES})
  openms_add_executable_compiler_flags(${_nuxl_test}_test)
  add_test(${_nuxl_test}_test ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_nuxl_test}_test)
  if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${_nuxl_test}_test PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endforeach()

# NuXLPresets.cpp parses the preset JSON; nlohmann::json is PRIVATE to libOpenMS, so it
# does not come with the OpenMS package either (see the tool's CMakeLists.txt).
if(NOT TARGET nlohmann_json::nlohmann_json)
  find_package(nlohmann_json CONFIG REQUIRED)
endif()
target_link_libraries(NuXLParameterParsing_test nlohmann_json::nlohmann_json)

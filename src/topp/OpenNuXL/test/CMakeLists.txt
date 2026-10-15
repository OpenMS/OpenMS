# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

# This directory is added by the OpenMS class-test project, which defines the class-test
# conventions in cmake/class_test_macros.cmake. Building it on its own is not supported:
# that needs OpenMS to export its test framework, which it deliberately does not (see
# src/testframework/CMakeLists.txt). State the dependency rather than leave a caller to
# discover it as an unknown command.
if(NOT COMMAND openms_add_tool_class_test)
  message(FATAL_ERROR
    "${CMAKE_CURRENT_SOURCE_DIR} is added by src/tests/class_tests/openms, which provides "
    "openms_add_tool_class_test(). It cannot be configured on its own.")
endif()

# Each entry lists the tool-local sources its test has to link: the class it covers plus the
# NuXL classes that class calls into.
set(_nuxl_NuXLDeisotoper_sources NuXLDeisotoper.cpp)
set(_nuxl_NuXLFragmentAdductDefinition_sources NuXLFragmentAdductDefinition.cpp)
set(_nuxl_NuXLModificationsGenerator_sources NuXLModificationsGenerator.cpp)
set(_nuxl_NuXLParameterParsing_sources NuXLParameterParsing.cpp NuXLFragmentAdductDefinition.cpp
                                       NuXLModificationsGenerator.cpp NuXLPresets.cpp)
set(_nuxl_NuXLReport_sources NuXLReport.cpp NuXLMarkerIonExtractor.cpp NuXLModificationsGenerator.cpp)

foreach(_nuxl_test NuXLDeisotoper NuXLFragmentAdductDefinition NuXLModificationsGenerator NuXLParameterParsing NuXLReport)
  set(_nuxl_test_sources)
  foreach(_nuxl_source ${_nuxl_${_nuxl_test}_sources})
    list(APPEND _nuxl_test_sources ${CMAKE_CURRENT_SOURCE_DIR}/../source/${_nuxl_source})
  endforeach()
  openms_add_tool_class_test(${_nuxl_test}_test
    SOURCES ${_nuxl_test_sources}
    INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)
endforeach()

# NuXLPresets.cpp parses the preset JSON; nlohmann::json is PRIVATE to libOpenMS, so it
# does not come with the OpenMS package either (see the tool's CMakeLists.txt).
if(NOT TARGET nlohmann_json::nlohmann_json)
  find_package(nlohmann_json CONFIG REQUIRED)
endif()
target_link_libraries(NuXLParameterParsing_test nlohmann_json::nlohmann_json)

# NuXLReport.cpp uses Boost.Math, and openms_add_tool_class_test() links ${OpenMS_LIBRARIES},
# which no longer carries Boost now that it is private to libOpenMS.
if(NOT TARGET Boost::boost)
  find_package(Boost 1.81 REQUIRED)
endif()
target_link_libraries(NuXLReport_test Boost::boost)

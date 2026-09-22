# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

# Class tests of the NuXL classes that live with this tool. The class-test project adds this
# directory and provides openms_add_tool_class_test(); see the block that does so in
# src/tests/class_tests/openms/CMakeLists.txt.
#
# Each entry lists the tool-local sources its test has to link: the class it covers plus the
# NuXL classes that class calls into.
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

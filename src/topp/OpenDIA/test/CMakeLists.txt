# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

if(COMMAND openms_add_tool_class_test)
  openms_add_tool_class_test(OpenDIAEntry_test)
else()
  # The boundary test needs only the public CLI SDK, not the class-test framework.
  add_executable(OpenDIAEntry_test OpenDIAEntry_test.cpp)
  add_test(NAME OpenDIAEntry_test COMMAND OpenDIAEntry_test)
endif()
target_link_libraries(OpenDIAEntry_test OpenDIA_impl)

if(COMMAND openms_add_tool_class_test)
  openms_add_tool_class_test(OpenDIARecords_test)
  target_compile_definitions(OpenDIARecords_test PRIVATE OPENDIA_CLASS_TEST)
else()
  add_executable(OpenDIARecords_test OpenDIARecords_test.cpp)
  add_test(NAME OpenDIARecords_test COMMAND OpenDIARecords_test)
endif()
target_link_libraries(OpenDIARecords_test OpenDIA_impl)

if(COMMAND openms_add_tool_class_test)
  openms_add_tool_class_test(OpenDIARecordsInterop_test)
else()
  add_executable(OpenDIARecordsInterop_test OpenDIARecordsInterop_test.cpp)
  add_test(NAME OpenDIARecordsInterop_test COMMAND OpenDIARecordsInterop_test)
endif()
# Use the real source-controlled OpenMS loader fixtures in both build modes.
get_filename_component(OPENDIA_TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../tests" ABSOLUTE)
target_compile_definitions(OpenDIARecordsInterop_test PRIVATE OPENDIA_TEST_DATA_DIR="${OPENDIA_TEST_DATA_DIR}")
target_link_libraries(OpenDIARecordsInterop_test OpenDIA_impl)

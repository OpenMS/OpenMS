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
  openms_add_tool_class_test(OpenDIAQuery_test)
  target_compile_definitions(OpenDIAQuery_test PRIVATE OPENDIA_CLASS_TEST)
else()
  add_executable(OpenDIAQuery_test OpenDIAQuery_test.cpp)
  add_test(NAME OpenDIAQuery_test COMMAND OpenDIAQuery_test)
endif()
target_link_libraries(OpenDIAQuery_test OpenDIA_impl)

if(COMMAND openms_add_tool_class_test)
  openms_add_tool_class_test(OpenDIAFamilies_test)
  target_compile_definitions(OpenDIAFamilies_test PRIVATE OPENDIA_CLASS_TEST)
else()
  add_executable(OpenDIAFamilies_test OpenDIAFamilies_test.cpp)
  add_test(NAME OpenDIAFamilies_test COMMAND OpenDIAFamilies_test)
endif()
target_link_libraries(OpenDIAFamilies_test OpenDIA_impl)

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

foreach(resource_test OpenDIAResourceBudget_test OpenDIABoundedQueue_test OpenDIAResourceInterop_test OpenDIARawAccess_test OpenDIARawFailure_test)
  if(COMMAND openms_add_tool_class_test)
    openms_add_tool_class_test(${resource_test})
    if(resource_test STREQUAL "OpenDIAResourceBudget_test" OR resource_test STREQUAL "OpenDIABoundedQueue_test")
      target_compile_definitions(${resource_test} PRIVATE OPENDIA_CLASS_TEST)
    endif()
  else()
    add_executable(${resource_test} ${resource_test}.cpp)
    add_test(NAME ${resource_test} COMMAND ${resource_test})
  endif()
  target_link_libraries(${resource_test} OpenDIA_impl)
  # A deadlock is a failing test, never an indefinitely stalled build.
  set_tests_properties(${resource_test} PROPERTIES TIMEOUT 30)
endforeach()
target_compile_definitions(OpenDIAResourceInterop_test PRIVATE OPENDIA_TEST_DATA_DIR="${OPENDIA_TEST_DATA_DIR}")
target_compile_definitions(OpenDIARawAccess_test PRIVATE OPENDIA_TEST_DATA_DIR="${OPENDIA_TEST_DATA_DIR}")
target_compile_definitions(OpenDIARawFailure_test PRIVATE OPENDIA_TEST_DATA_DIR="${OPENDIA_TEST_DATA_DIR}")

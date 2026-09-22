# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

# Class tests of the workflow-specific classes that live with this tool. The class-test
# project adds this directory and provides openms_add_tool_class_test(); see the OpenNuXL
# block in src/tests/class_tests/openms/CMakeLists.txt. Each class is standalone, so each
# test links only the class it covers.
foreach(_ms1_test MS1LabeledFAIMS MS1LabeledRatioQuantifier MS1LabeledSpectra)
  openms_add_tool_class_test(${_ms1_test}_test
    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../${_ms1_test}.cpp
    INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..)
endforeach()

# Optional Thermo RAW fixture for MS1LabeledSpectra::load(); the test skips without it.
if (WITH_THERMO_RAW AND THERMO_RAW_TEST_DATA)
  file(TO_CMAKE_PATH "${THERMO_RAW_TEST_DATA}" _ms1_thermo_test_data)
  target_compile_definitions(MS1LabeledSpectra_test PRIVATE
    THERMO_RAW_TEST_DATA="${_ms1_thermo_test_data}")
  set_tests_properties(MS1LabeledSpectra_test PROPERTIES TIMEOUT 1200)
  if (COMMAND openms_thermo_bridge_copy_runtime_files)
    openms_thermo_bridge_copy_runtime_files(TARGET MS1LabeledSpectra_test)
  endif()
endif()

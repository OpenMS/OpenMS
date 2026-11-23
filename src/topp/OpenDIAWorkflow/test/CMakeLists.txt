# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

add_executable(OpenDIASearch_test OpenDIASearch_test.cpp)
target_link_libraries(OpenDIASearch_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIASearch_test COMMAND OpenDIASearch_test)

add_executable(OpenDIATraces_test OpenDIATraces_test.cpp)
target_link_libraries(OpenDIATraces_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIATraces_test COMMAND OpenDIATraces_test ${CMAKE_CURRENT_BINARY_DIR}/OpenDIATraces_test)

add_executable(OpenDIACalibration_test OpenDIACalibration_test.cpp)
target_link_libraries(OpenDIACalibration_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIACalibration_test COMMAND OpenDIACalibration_test)

add_executable(OpenDIAStage2_test OpenDIAStage2_test.cpp)
target_link_libraries(OpenDIAStage2_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAStage2_test COMMAND OpenDIAStage2_test)

add_executable(OpenDIAPeakPicker_test OpenDIAPeakPicker_test.cpp)
target_link_libraries(OpenDIAPeakPicker_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAPeakPicker_test COMMAND OpenDIAPeakPicker_test)

add_executable(OpenDIAInterference_test OpenDIAInterference_test.cpp)
target_link_libraries(OpenDIAInterference_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAInterference_test COMMAND OpenDIAInterference_test)

add_executable(OpenDIAScorer_test OpenDIAScorer_test.cpp)
target_link_libraries(OpenDIAScorer_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAScorer_test COMMAND OpenDIAScorer_test)

add_executable(OpenDIAGBM_test OpenDIAGBM_test.cpp)
target_link_libraries(OpenDIAGBM_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAGBM_test COMMAND OpenDIAGBM_test)

add_executable(OpenDIAFragmentGrid_test OpenDIAFragmentGrid_test.cpp)
target_link_libraries(OpenDIAFragmentGrid_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAFragmentGrid_test COMMAND OpenDIAFragmentGrid_test)
add_executable(OpenDIAMzIndex_test OpenDIAMzIndex_test.cpp)
target_link_libraries(OpenDIAMzIndex_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIAMzIndex_test COMMAND OpenDIAMzIndex_test)
add_executable(OpenDIARTRefine_test OpenDIARTRefine_test.cpp)
target_link_libraries(OpenDIARTRefine_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIARTRefine_test COMMAND OpenDIARTRefine_test)
# Design rule: no OpenSWATH code in the OpenDIA core library.
add_test(NAME OpenDIA_no_openswath
  COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/.. -P ${CMAKE_CURRENT_SOURCE_DIR}/OpenDIA_no_openswath.cmake)

if(EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML" AND
   EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp")
  # OpenDIAWorkflow fails closed: a library too small to calibrate, and two inputs with one file name (outputs
  # would collide).
  add_test(NAME OpenDIAWorkflow_fails_closed_on_tiny_library
    COMMAND $<TARGET_FILE:OpenDIAWorkflow>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -out_dir ${CMAKE_CURRENT_BINARY_DIR}/OpenDIAWorkflow_tiny)
  set_tests_properties(OpenDIAWorkflow_fails_closed_on_tiny_library PROPERTIES WILL_FAIL TRUE)
  add_test(NAME OpenDIAWorkflow_rejects_duplicate_run_names
    COMMAND $<TARGET_FILE:OpenDIAWorkflow>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -out_dir ${CMAKE_CURRENT_BINARY_DIR}/OpenDIAWorkflow_duplicate)
  set_tests_properties(OpenDIAWorkflow_rejects_duplicate_run_names PROPERTIES WILL_FAIL TRUE)
endif()

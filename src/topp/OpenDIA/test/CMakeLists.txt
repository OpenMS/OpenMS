# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

add_executable(OpenDIASearch_test OpenDIASearch_test.cpp)
target_link_libraries(OpenDIASearch_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIASearch_test COMMAND OpenDIASearch_test)

add_executable(OpenDIARerank_test OpenDIARerank_test.cpp)
target_link_libraries(OpenDIARerank_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIARerank_test COMMAND OpenDIARerank_test)

add_executable(OpenDIATraces_test OpenDIATraces_test.cpp)
target_link_libraries(OpenDIATraces_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIATraces_test COMMAND OpenDIATraces_test ${CMAKE_CURRENT_BINARY_DIR}/OpenDIATraces_test)

add_executable(OpenDIACalibration_test OpenDIACalibration_test.cpp)
target_link_libraries(OpenDIACalibration_test PRIVATE OpenDIA_impl)
add_test(NAME OpenDIACalibration_test COMMAND OpenDIACalibration_test)

if(TARGET OpenDIA_cnn)
  add_executable(OpenDIACNN_test OpenDIACNN_test.cpp)
  target_link_libraries(OpenDIACNN_test PRIVATE OpenDIA_cnn)
  add_test(NAME OpenDIACNN_test COMMAND OpenDIACNN_test)
  add_executable(OpenDIACNNQValues_test OpenDIACNNQValues_test.cpp)
  target_link_libraries(OpenDIACNNQValues_test PRIVATE OpenDIA_cnn)
  add_test(NAME OpenDIACNNQValues_test COMMAND OpenDIACNNQValues_test)
  add_executable(OpenDIACNNRTWindow_test OpenDIACNNRTWindow_test.cpp)
  target_link_libraries(OpenDIACNNRTWindow_test PRIVATE OpenDIA_cnn)
  add_test(NAME OpenDIACNNRTWindow_test COMMAND OpenDIACNNRTWindow_test)
endif()

if(EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML" AND
   EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp")
  if(TARGET OpenDIACNN)
    add_test(NAME OpenDIACNN_cli_safety
      COMMAND ${CMAKE_COMMAND} -DBINARY=$<TARGET_FILE:OpenDIACNN>
        -DLIBRARY=${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/OpenDIACNN_cli_safety
        -P ${CMAKE_CURRENT_SOURCE_DIR}/OpenDIACNN_cli_test.cmake)
  endif()
  # OpenDIACalibrate refuses to overwrite its inputs, and fails closed (no silent fallback) when the tiny
  # fixture library has too few target-decoy pairs to calibrate.
  add_test(NAME OpenDIACalibrate_refuses_input_overwrite
    COMMAND $<TARGET_FILE:OpenDIACalibrate>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -out ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp)
  set_tests_properties(OpenDIACalibrate_refuses_input_overwrite PROPERTIES WILL_FAIL TRUE)
  add_test(NAME OpenDIACalibrate_fails_closed_on_tiny_library
    COMMAND $<TARGET_FILE:OpenDIACalibrate>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -out ${CMAKE_CURRENT_BINARY_DIR}/OpenDIACalibrate_tiny.json)
  set_tests_properties(OpenDIACalibrate_fails_closed_on_tiny_library PROPERTIES WILL_FAIL TRUE)
  # Panel escalation stays fail-closed when the library cannot supply a larger panel, and rejects a cap below the start.
  add_test(NAME OpenDIACalibrate_escalation_fails_closed_on_tiny_library
    COMMAND $<TARGET_FILE:OpenDIACalibrate>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -max_panel_targets 40000
      -out ${CMAKE_CURRENT_BINARY_DIR}/OpenDIACalibrate_tiny_escalation.json)
  set_tests_properties(OpenDIACalibrate_escalation_fails_closed_on_tiny_library PROPERTIES WILL_FAIL TRUE)
  add_test(NAME OpenDIACalibrate_rejects_small_escalation_cap
    COMMAND $<TARGET_FILE:OpenDIACalibrate>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -max_panel_targets 100
      -out ${CMAKE_CURRENT_BINARY_DIR}/OpenDIACalibrate_small_cap.json)
  set_tests_properties(OpenDIACalibrate_rejects_small_escalation_cap PROPERTIES WILL_FAIL TRUE)
  # OpenDIA -calibration rejects an invalid calibration file instead of silently extracting uncalibrated.
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/OpenDIA_invalid_calibration.json "{\"version\": 1, \"rt\": {}}\n")
  add_test(NAME OpenDIA_rejects_invalid_calibration
    COMMAND $<TARGET_FILE:OpenDIA>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -calibration ${CMAKE_CURRENT_BINARY_DIR}/OpenDIA_invalid_calibration.json
      -out ${CMAKE_CURRENT_BINARY_DIR}/OpenDIA_invalid_calibration.tsv)
  set_tests_properties(OpenDIA_rejects_invalid_calibration PROPERTIES WILL_FAIL TRUE)
  add_test(NAME OpenDIA_streaming_smoke
    COMMAND $<TARGET_FILE:OpenDIA>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -out ${CMAKE_CURRENT_BINARY_DIR}/OpenDIA_streaming_smoke.tsv
      -min_transitions 1)
  add_test(NAME OpenDIA_feature_table_smoke
    COMMAND $<TARGET_FILE:OpenDIA>
      -in ${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML
      -lib ${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp
      -out ${CMAKE_CURRENT_BINARY_DIR}/OpenDIA_feature_table_smoke.tsv
      -features_out ${CMAKE_CURRENT_BINARY_DIR}/OpenDIA_feature_table_smoke.features.tsv
      -min_transitions 1)
endif()

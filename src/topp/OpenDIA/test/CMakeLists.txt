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

if(EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/OpenSwathWorkflow_17_input.mzML" AND
   EXISTS "${CMAKE_SOURCE_DIR}/src/tests/topp/TargetedFileConverter_10_input.pqp")
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

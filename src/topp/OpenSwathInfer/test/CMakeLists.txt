# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

if(NOT COMMAND openms_add_tool_class_test)
  message(FATAL_ERROR
    "${CMAKE_CURRENT_SOURCE_DIR} is added by src/tests/class_tests/openms, which provides "
    "openms_add_tool_class_test(). It cannot be configured on its own.")
endif()

openms_add_tool_class_test(OpenSwathInference_test
  SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathGeneInference.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathPeptideInference.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathProteinInference.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

# PeptidoformInference_test was gated on OpenSWATH before the move.
if(NOT DISABLE_OPENSWATH)
  openms_add_tool_class_test(PeptidoformInference_test
    SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/../source/OpenSwathPeptidoformInference.cpp
    INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../include)

endif()

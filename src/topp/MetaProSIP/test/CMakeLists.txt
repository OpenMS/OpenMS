# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

# This directory is added by the OpenMS class-test project (see src/topp/MS1LabeledWorkflow/test).
if(NOT COMMAND openms_add_tool_class_test)
  message(FATAL_ERROR
    "${CMAKE_CURRENT_SOURCE_DIR} is added by src/tests/class_tests/openms, which provides "
    "openms_add_tool_class_test(). It cannot be configured on its own.")
endif()

openms_add_tool_class_test(MetaProSIPPseudoFeatures_test
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../MetaProSIPPseudoFeatures.cpp
  INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..)

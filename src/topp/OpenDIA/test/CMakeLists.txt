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

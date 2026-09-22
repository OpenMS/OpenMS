# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

# The class-test conventions, for classes that live with a TOPP tool rather than in the
# library. A tool's test/ folder calls openms_add_tool_class_test() and states only what is
# specific to it: which sources its test links and where the tool's headers are.
#
# This lives here rather than inside src/tests/class_tests/openms/CMakeLists.txt so that a
# tool's test/ folder has one documented entry point instead of depending on the internals
# of the project that happens to add it, and so the same entry point can be installed with
# the OpenMS package once the test framework is exported. It is not installed today: the
# package has no TestFramework component, so there is nothing out of tree to link against
# (see src/testframework/CMakeLists.txt for why the target is not exported).
#
# The caller sets OPENMS_TEST_SUPPORT_SOURCE. OpenMSTestSupport.cpp registers libOpenMS
# behavior (unique-ID seeding, exception naming) with the standard-library-only test
# framework, and every class test links it; it belongs to the class-test project today, so
# an out-of-tree consumer has nothing to point this at yet.

# openms_add_tool_class_test(<name>_test
#                            SOURCES <tool sources the test links>
#                            INCLUDE_DIRS <dirs holding the tool's headers>)
#
# Builds <name>_test from <name>_test.cpp in the calling directory and registers it, the way
# the loop in src/tests/class_tests/openms/CMakeLists.txt builds a library class test.
function(openms_add_tool_class_test _name)
  cmake_parse_arguments(_tct "" "" "SOURCES;INCLUDE_DIRS" ${ARGN})
  if(_tct_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "openms_add_tool_class_test(${_name}): unexpected arguments "
                        "'${_tct_UNPARSED_ARGUMENTS}'")
  endif()
  if(NOT OPENMS_TEST_SUPPORT_SOURCE)
    message(FATAL_ERROR
      "openms_add_tool_class_test(${_name}): OPENMS_TEST_SUPPORT_SOURCE is not set. It has to "
      "name OpenMSTestSupport.cpp, which every class test links; without it, reference-file "
      "comparisons of ID-bearing output become nondeterministic.")
  endif()

  add_executable(${_name} ${_name}.cpp ${_tct_SOURCES} "${OPENMS_TEST_SUPPORT_SOURCE}")
  target_include_directories(${_name} PRIVATE ${_tct_INCLUDE_DIRS})
  target_link_libraries(${_name} OpenMSTestFramework ${OpenMS_LIBRARIES})
  if(COMMAND openms_add_executable_compiler_flags)
    openms_add_executable_compiler_flags(${_name})
  endif()
  add_test(${_name} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_name})
  # only add OPENMP flags to gcc linker (except Mac OS X, due to a compiler bug)
  if(OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${_name} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endfunction()

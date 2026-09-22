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
# of the project that happens to add it.
#
# In-tree only, and not installed. OpenMS deliberately does not export its test framework
# (see src/testframework/CMakeLists.txt), so a tool's test/ folder cannot be configured
# outside this build; those folders therefore assert on this command rather than offer a
# fallback that could not succeed. Making them buildable out of tree needs the framework
# exported as a package component -- a decision in its own right, which this file does not
# anticipate.
#
# The OpenMSTestSupport target registers libOpenMS behavior (unique-ID seeding, exception
# naming) with the standard-library-only test framework. Every class test links it, and
# linking it is the only thing a test has to do to get it: see src/testframework.

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

  add_executable(${_name} ${_name}.cpp ${_tct_SOURCES})
  target_include_directories(${_name} PRIVATE ${_tct_INCLUDE_DIRS})
  target_link_libraries(${_name} OpenMSTestSupport ${OpenMS_LIBRARIES})
  if(COMMAND openms_add_executable_compiler_flags)
    openms_add_executable_compiler_flags(${_name})
  endif()
  add_test(${_name} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_name})
  # only add OPENMP flags to gcc linker (except Mac OS X, due to a compiler bug)
  if(OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${_name} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endfunction()

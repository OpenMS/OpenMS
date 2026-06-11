# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


# CMake sub-project for OpenMS code examples

project("OpenMS_code_examples")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(EXAMPLES_executables)

# ensure the tutorials go into bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# get the actual code examples
include(${CMAKE_CURRENT_LIST_DIR}/executables.cmake)


# add the targets
foreach(i ${EXAMPLES_executables})
  add_executable(${i} ${CMAKE_CURRENT_LIST_DIR}/${i}.cpp)
  target_link_libraries(${i} OpenMS)

  # if not GUI and no style tests requested -> add as test
  if(NOT ENABLE_STYLE_TESTING)
    add_test(NAME ${i}
             COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i})
  endif()
endforeach(i)
## expect "Tutorial_Template"" to fail, since its a proper TOPP tool invoked without arguments
## i.e. TOPPBase will print a usage() and return an error code
if(NOT ENABLE_STYLE_TESTING)
  set_tests_properties(Tutorial_Template PROPERTIES WILL_FAIL TRUE)
endif()

# add collection target
add_custom_target(Tutorials_build)
add_dependencies(Tutorials_build ${EXAMPLES_executables})

if(WITH_GUI)
  # add the GUI tutorials -- no tests
  foreach(i ${GUI_EXAMPLES_executables})
    add_executable(${i} ${CMAKE_CURRENT_LIST_DIR}/${i}.cpp)
    target_link_libraries(${i} OpenMS OpenMS_GUI)
  endforeach(i)

  add_dependencies(Tutorials_build ${GUI_EXAMPLES_executables})
endif()

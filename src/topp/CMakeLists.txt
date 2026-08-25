# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS_TOPP")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# --------------------------------------------------------------------------
# OpenMS' TOPP tools

# add all the tools. We support two modes
# 1. either as single MyToolName.cpp file
# 2. via a subfolder MyToolName/ with MyToolName/MyToolName.cpp and (optionally) additional header/source files
#    Use this mode for classes which are so specific to a single tool that they do not belong
#    into the OpenMS library (and e.g. should not be wrapped in pyOpenMS).
include(executables.cmake)

#------------------------------------------------------------------------------
# Support for class tests living in a tool's subdirectory (tests/): a shared
# OpenMS/test_config.h (whose data path points to the central class test data
# directory) plus a macro to register such a test with CTest.
set(CF_OPENMS_TEST_DATA_PATH "${OpenMS_TOPP_SOURCE_DIR}/../tests/class_tests/openms/data/")
configure_file(${OpenMS_TOPP_SOURCE_DIR}/../tests/class_tests/openms/include/OpenMS/test_config.h.in
               ${OpenMS_TOPP_BINARY_DIR}/tool_tests_include/OpenMS/test_config.h)

# Registers tests/<test_name>.cpp of the calling tool subdirectory as a class
# test executable linked against the tool's helper library <tool_lib>. The
# tests honor the same switches as the library class tests. Note:
# ENABLE_CLASS_TESTING is declared in src/tests, which is processed after
# src/topp - hence the 'NOT DEFINED' fallback to its default ON.
macro(openms_add_tool_class_test test_name tool_lib)
  if((NOT DEFINED ENABLE_CLASS_TESTING OR ENABLE_CLASS_TESTING) AND "${PACKAGE_TYPE}" STREQUAL "none")
    add_executable(${test_name} tests/${test_name}.cpp)
    openms_add_executable_compiler_flags(${test_name})
    target_link_libraries(${test_name} ${tool_lib})
    target_include_directories(${test_name} PRIVATE ${OpenMS_TOPP_BINARY_DIR}/tool_tests_include)
    # keep the test binaries out of the bin/ directory holding the TOPP tools
    set_target_properties(${test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
    if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set_target_properties(${test_name} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
    endif()
  endif()
endmacro()

foreach(i ${TOPP_executables})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${i})
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${i})
#[[
# CMakeLists.txt in the subdirectory for MyToolName could look like this:
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(NameOfTheTool)

# allow building separately from OpenMS
if(NOT TARGET OpenMS)
  find_package(OpenMS)
endif()

set(tool_sources
  additional_file1.cpp
  additional_file2.cpp
  additional_header1.h
  additional_header2.h
)
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${tool_sources}) # Add the main cpp file and additional files
openms_add_executable_compiler_flags(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} OpenMS) # Link the executable with OpenMS

# If the additional classes have class tests, compile the classes into a static helper
# library '${PROJECT_NAME}_lib' instead, link both the tool and the tests against it, and
# put the tests into a 'tests/' subfolder next to the tool: since src/topp is built on top
# of libOpenMS, the tool registers its own class tests with CTest instead of adding them
# to src/tests/class_tests/openms (see src/topp/DatabaseSuitability or src/topp/CometAdapter
# for the full pattern, including the ENABLE_CLASS_TESTING/PACKAGE_TYPE guard):
add_library(${PROJECT_NAME}_lib STATIC additional_file1.cpp ...)
openms_add_library_compiler_flags(${PROJECT_NAME}_lib)
target_include_directories(${PROJECT_NAME}_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME}_lib PUBLIC OpenMS)
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
openms_add_executable_compiler_flags(${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_lib)
openms_add_tool_class_test(MyClass_test ${PROJECT_NAME}_lib) # expects tests/MyClass_test.cpp
#]]
  else()
    add_executable(${i} ${i}.cpp)
    target_link_libraries(${i} OpenMS)
    openms_add_executable_compiler_flags(${i})
  endif()  
  
  # we also want to install each topp tool
  install_tool(${i})
  if (OPENMP_FOUND AND NOT MSVC)
    set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endforeach(i)

# special requirements of some tools
target_link_libraries(FileMerger Boost::regex)
target_link_libraries(IsobaricWorkflow Eigen3::Eigen)
target_link_libraries(MetaProSIP Eigen3::Eigen)


# collection target
add_custom_target(TOPP)
add_dependencies(TOPP ${TOPP_executables})

## export the list of TOPP tools into CACHE
set(TOPP_TOOLS ${TOPP_executables}
    CACHE INTERNAL "OpenMS' TOPP tools" FORCE)

openms_doc_path("${PROJECT_SOURCE_DIR}")

# Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
include(executables.cmake)
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
target_link_libraries(${PROJECT_NAME} OpenMS) # Link the executable with OpenMS
#]] 
  else()
    add_executable(${i} ${i}.cpp)
    target_link_libraries(${i} OpenMS)
  endif()  
  
  # we also want to install each topp tool
  install_tool(${i})
  if (OPENMP_FOUND AND NOT MSVC)
    set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endforeach(i)

# some regular TOPP tools need the GUI lib, only build them when WITH_GUI is enabled
if(WITH_GUI)
  foreach(i ${TOPP_executables_with_GUIlib})
    add_executable(${i} ${i}.cpp)
    # we also want to install each topp tool
    install_tool(${i})
    target_link_libraries(${i} OpenMS OpenMS_GUI)
  endforeach(i)

  set(TOPP_executables ${TOPP_executables} ${TOPP_executables_with_GUIlib})
endif()

# special requirements of some tools
target_link_libraries(FileMerger Boost::regex)
target_link_libraries(DecoyDatabase Boost::regex)
target_link_libraries(NucleicAcidSearchEngine Boost::regex)

# collection target
add_custom_target(TOPP)
add_dependencies(TOPP ${TOPP_executables})

## export the list of TOPP tools into CACHE
set(TOPP_TOOLS ${TOPP_executables}
    CACHE INTERNAL "OpenMS' TOPP tools" FORCE)

openms_doc_path("${PROJECT_SOURCE_DIR}")

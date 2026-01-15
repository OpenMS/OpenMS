# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_class_tests_openms_gui")

set(visual_executables_list
  AxisTickCalculator_test
  GUIHelpers_test
  MultiGradient_test
)

set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# TODO should not be required anymore
#------------------------------------------------------------------------------
# required to build GUI tests (e.g., TOPPView_test)
include_directories(SYSTEM ${OpenMS_GUI_INCLUDE_DIRECTORIES})

# add include to find moced files
include_directories(${PROJECT_BINARY_DIR})

# required for GUI tests to find the ui_foo.h files (e.g. TSGDialog_test)
include_directories(${CMAKE_BINARY_DIR}/src/openms_gui)

#------------------------------------------------------------------------------
# VISUAL
add_custom_target(VISUAL_TEST)
add_dependencies(VISUAL_TEST ${visual_executables_list})

foreach(_class_test ${visual_executables_list})
  add_executable(${_class_test} source/${_class_test}.cpp)
  target_link_libraries(${_class_test} ${OpenMS_GUI_LIBRARIES})
  openms_add_executable_compiler_flags(${_class_test})
  add_test(${_class_test} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_class_test})
endforeach(_class_test)

#------------------------------------------------------------------------------
# GUI tests special treatment - Apply MOC compiler to GUI
#------------------------------------------------------------------------------

# we only want to run those tests if we have a running XServer
if(HAS_XSERVER)
  set(GUI_executables_list
    TOPPView_test
    TSGDialog_test
  )

  #------------------------------------------------------------------------------
  # Use QtTest lib

  find_package(Qt6 REQUIRED
    Test
  )
  if (NOT Qt6Test_FOUND)
    message(STATUS "QtTest module not found!")
    message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
  endif()

  #------------------------------------------------------------------------------
  # Add GUI tests
  foreach(_gui_test ${GUI_executables_list})
    # add the executable (Note: The header needs to be added as source file, too. Otherwise automoc will not find it.
    add_executable(${_gui_test} source/GUI/${_gui_test}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/include/GUI/${_gui_test}.h)
    # add headers
    target_include_directories(${_gui_test} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/GUI)
    # link against openms_gui and qt
    target_link_libraries(${_gui_test} ${OpenMS_GUI_LIBRARIES} ${QT_LIBRARIES} Qt6::Test)
    # add the test
    add_test(${_gui_test} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_gui_test})
  endforeach(_gui_test)

  #------------------------------------------------------------------------------
  # GUI tests
  add_custom_target(GUI_TEST)
  add_dependencies(GUI_TEST ${GUI_executables_list})

endif()

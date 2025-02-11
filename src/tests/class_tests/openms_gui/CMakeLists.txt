# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright OpenMS Inc. -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-present.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

# CMake sub-project for OpenMS code examples

project("OpenMS_code_examples")
cmake_minimum_required(VERSION 3.8.0 FATAL_ERROR)

set(EXAMPLES_executables)

# ensure the tutorials got into bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# get the actual code examples
include(${CMAKE_CURRENT_LIST_DIR}/executables.cmake)

# add OpenMS includes
include_directories(SYSTEM ${OpenMS_INCLUDE_DIRECTORIES})

#------------------------------------------------------------------------------
# QT dependencies
find_package(Qt5 COMPONENTS Core Network Sql REQUIRED)

# add the targets
foreach(i ${EXAMPLES_executables})
  add_executable(${i} ${CMAKE_CURRENT_LIST_DIR}/${i}.cpp)
  target_link_libraries(${i} ${OpenMS_LIBRARIES})

  # if not GUI and no style tests requested -> add as test
  if(NOT ENABLE_STYLE_TESTING)
    add_test(NAME ${i}
             COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${i} ${CMAKE_CURRENT_LIST_DIR})
  endif()
endforeach(i)

# add collection target
add_custom_target(Tutorials_build)
add_dependencies(Tutorials_build ${EXAMPLES_executables})

if(WITH_GUI)
  find_package(Qt5 COMPONENTS Core Network OpenGL Svg REQUIRED)
  # add OpenMS GUI includes
  include_directories(SYSTEM ${OpenMS_GUI_INCLUDE_DIRECTORIES})

  # add the GUI tutorials -- no tests
  foreach(i ${GUI_EXAMPLES_executables})
    add_executable(${i} ${CMAKE_CURRENT_LIST_DIR}/${i}.cpp)
    target_link_libraries(${i} ${OpenMS_GUI_LIBRARIES})
  endforeach(i)

  add_dependencies(Tutorials_build ${GUI_EXAMPLES_executables})
endif()

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

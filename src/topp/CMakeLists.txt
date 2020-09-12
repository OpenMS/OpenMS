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
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS_TOPP")
cmake_minimum_required(VERSION 3.8.0 FATAL_ERROR)

# --------------------------------------------------------------------------
# OpenMS' TOPP tools

# add OpenMS includes
include_directories(SYSTEM ${OpenMS_INCLUDE_DIRECTORIES})
add_definitions(/DBOOST_ALL_NO_LIB)

find_package(Qt5 COMPONENTS Core Network Sql REQUIRED)

# add all the tools
include(executables.cmake)
foreach(i ${TOPP_executables})
  add_executable(${i} ${i}.cpp)
  target_link_libraries(${i} ${OpenMS_LIBRARIES})
  # we also want to install each topp tool
  install_tool(${i})
  if (OPENMP_FOUND AND NOT MSVC)
    set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endforeach(i)

# some regular TOPP tools need the GUI lib, only build them when WITH_GUI is enabled
if(WITH_GUI)
  ## TODO Check if the find command actually does anything but double check.
  ## Libraries are I think anyway already hardcoded in OpenMS_GUI_LIBRARIES
  find_package(Qt5 COMPONENTS Core Network Widgets Svg OpenGL PrintSupport Sql REQUIRED)
  if (NOT Qt5Widgets_FOUND)
    message(STATUS "QtWidgets module not found!")
    message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
  endif()

  ## some regular TOPP tools need the GUI lib
  include_directories(SYSTEM ${OpenMS_GUI_INCLUDE_DIRECTORIES})
  foreach(i ${TOPP_executables_with_GUIlib})
    add_executable(${i} ${i}.cpp)
    # we also want to install each topp tool
    install_tool(${i})
    target_link_libraries(${i} ${OpenMS_GUI_LIBRARIES})
  endforeach(i)

  set(TOPP_executables ${TOPP_executables} ${TOPP_executables_with_GUIlib})
endif()

# collection target
add_custom_target(TOPP)
add_dependencies(TOPP ${TOPP_executables})

## export the list of TOPP tools into CACHE
set(TOPP_TOOLS ${TOPP_executables}
    CACHE INTERNAL "OpenMS' TOPP tools" FORCE)

openms_doc_path("${PROJECT_SOURCE_DIR}")

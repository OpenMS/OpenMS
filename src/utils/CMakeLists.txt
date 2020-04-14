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

project("OpenMS_UTILS")
cmake_minimum_required(VERSION 3.8.0 FATAL_ERROR)

# --------------------------------------------------------------------------
# OpenMS' UTILS tools

# add OpenMS includes
include_directories(SYSTEM ${OpenMS_INCLUDE_DIRECTORIES})
add_definitions(/DBOOST_ALL_NO_LIB)

find_package(Qt5 COMPONENTS Core Network Sql REQUIRED)

# add all the tools
set(UTILS_executables)
include(executables.cmake)
foreach(i ${UTILS_executables})
  add_executable(${i} ${i}.cpp)
  # we also want to install each util
  install_tool(${i})
  target_link_libraries(${i} ${OpenMS_LIBRARIES})
endforeach(i)

# some regular UTILS tools need the GUI lib, only build them when WITH_GUI is enabled
if(WITH_GUI)
  find_package(Qt5 COMPONENTS Core Network Widgets Svg OpenGL Sql REQUIRED)
  if (NOT Qt5Widgets_FOUND)
    message(STATUS "QtWidgets module not found!")
    message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
  endif()
  ## some regular UTILS need the GUI lib
  include_directories(SYSTEM ${OpenMS_GUI_INCLUDE_DIRECTORIES})
  foreach(i ${UTILS_executables_with_GUIlib})
    add_executable(${i} ${i}.cpp)
    # we also want to install each util
    install_tool(${i})
    target_link_libraries(${i} ${OpenMS_GUI_LIBRARIES})
  endforeach()

  set(UTILS_executables ${UTILS_executables} ${UTILS_executables_with_GUIlib})
endif()

# collection target
add_custom_target(UTILS)
add_dependencies(UTILS ${UTILS_executables})

# FFSuperHirn requires linking against SuperHirn and pthread (needs find_package (Threads REQUIRED))
include_directories(SYSTEM ${SuperHirn_INCLUDE_DIRECTORIES})
target_link_libraries(FeatureFinderSuperHirn SuperHirn ${CMAKE_THREAD_LIBS_INIT})

## export the list of UTILS tools into CACHE
set(UTILS_TOOLS ${UTILS_executables}
    CACHE INTERNAL "OpenMS' UTILS tools" FORCE)

openms_doc_path("${PROJECT_SOURCE_DIR}")

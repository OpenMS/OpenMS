# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

project("OpenMS_GUI")
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

# --------------------------------------------------------------------------
# Find Qt
#------------------------------------------------------------------------------

find_package(Qt4 REQUIRED QtGui QtOpenGL QtSVG QtWebKit QtNetwork)
if (NOT QT4_FOUND)
  message(STATUS "Qt4 not found!")
  message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
endif()
# required to indicate we use qtopengl
set(QT_USE_QTOPENGL TRUE)
include(${QT_USE_FILE})

# --------------------------------------------------------------------------
# include custom macros for openms_gui lib
include(qt_wrap_ui.cmake)
include(add_mac_bundle.cmake)

#------------------------------------------------------------------------------
# Find boost
find_boost()
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost was not found!")
endif()

# --------------------------------------------------------------------------
# big include file for headers and cpp files
include (${PROJECT_SOURCE_DIR}/includes.cmake)

# --------------------------------------------------------------------------
# collect link dependencies
set(OpenMS_GUI_DEP_LIBRARIES OpenMS
                             ${QT_LIBRARIES})

# --------------------------------------------------------------------------
# handle OpenGL on macosx systems
if(APPLE)
  find_package(OpenGL REQUIRED)
  set(OpenMS_GUI_DEP_LIBRARIES ${OpenMS_GUI_DEP_LIBRARIES}
                               ${OPENGL_gl_LIBRARY}
                               ${OPENGL_glu_LIBRARY})
endif()

#------------------------------------------------------------------------------
# add the library
openms_add_library(TARGET_NAME OpenMS_GUI
                   SOURCE_FILES ${OpenMSVisual_sources}
                   HEADER_FILES ${OpenMSVisual_sources_h}
                   INTERNAL_INCLUDES ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include
                   EXTERNAL_INCLUDES ${OpenMS_INCLUDE_DIRECTORIES}
                   LINK_LIBRARIES ${OpenMS_GUI_DEP_LIBRARIES}
                   DLL_EXPORT_PATH "OpenMS/VISUAL/")

# --------------------------------------------------------------------------
# additional linker flags required by openms_gui
if (MSVC)
  set (GUI_lnk_flags "/FORCE:MULTIPLE /INCREMENTAL:NO /ignore:4006 /ignore:4088")
  set_target_properties(OpenMS_GUI PROPERTIES LINK_FLAGS ${GUI_lnk_flags}) ## allow multiple definitions of symbols (e.g. from template instanciations or STL-derived classes)
endif()

#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
openms_doc_path("${PROJECT_SOURCE_DIR}/include")
openms_doc_path("${PROJECT_SOURCE_DIR}/source/VISUAL/APPLICATIONS")

# --------------------------------------------------------------------------
# GUI tools
set(GUI_DIR source/VISUAL/APPLICATIONS/GUITOOLS)
include(${GUI_DIR}/executables.cmake)
foreach(i ${GUI_executables})
  set(resource_file ${PROJECT_SOURCE_DIR}/${GUI_DIR}/${i}.rc)
  set(resource_dir ${PROJECT_SOURCE_DIR}/${GUI_DIR}/${i}-resources/)
  ## add icons to TOPPView and INIFileEditor
  if (MSVC AND EXISTS ${resource_file})
    message(STATUS "Setting resource file ${resource_file} for ${i}")
    add_executable(${i} ${GUI_DIR}/${i}.cpp ${resource_file})
    install_tool(${i})
  elseif(APPLE AND EXISTS ${resource_dir})
    add_mac_app_bundle(${i})
  else()
    if(APPLE)
      message(STATUS "No icon file (${i}.icns) found in ${resource_dir}. Will not build ${i} as app bundle.")
    else()
      message(STATUS "No resource file (${resource_file}) found for ${i}. No icon will be embedded.")
    endif()
    add_executable(${i} ${GUI_DIR}/${i}.cpp)
    install_tool(${i})
  endif()

  ## append visual lib as dependency for GUI tools
  target_link_libraries(${i} ${OpenMS_GUI_LIBRARIES})
  if (OPENMP_FOUND AND NOT MSVC)
    set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endforeach(i)
add_custom_target(GUI)
add_dependencies(GUI ${GUI_executables})

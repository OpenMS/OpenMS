# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-.
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
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(CMAKE_AUTOUIC OFF) # leave off since this is broken, see qt5_wrap_ui below
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# --------------------------------------------------------------------------
# include custom macros for openms_gui lib
#include(qt_wrap_ui.cmake)
include(add_mac_bundle.cmake)

# --------------------------------------------------------------------------
# big include file for headers and cpp files
include (${PROJECT_SOURCE_DIR}/includes.cmake)

# --------------------------------------------------------------------------

# manually extract .ui files and track them (trigger rebuild of ui_CLASS.h upon changing CLASS.ui)
# see https://gitlab.kitware.com/cmake/cmake/issues/18741
set(UI_SOURCES ${OpenMSVisual_sources}) # copy list of .cpp + .ui files
list(FILTER UI_SOURCES INCLUDE REGEX "\\.ui$")
qt5_wrap_ui(UI_SOURCES2 ${UI_SOURCES}) # transform .ui to ui_....h files and add them to OpenMSVisual_sources_h

#------------------------------------------------------------------------------
# add the library
openms_add_library(TARGET_NAME OpenMS_GUI
                   SOURCE_FILES ${OpenMSVisual_sources}
                   HEADER_FILES ${OpenMSVisual_sources_h}
                   INTERNAL_INCLUDES ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include
                   EXTERNAL_INCLUDES
                   LINK_LIBRARIES ${OpenMS_GUI_DEP_LIBRARIES}
                   DLL_EXPORT_PATH "OpenMS/VISUAL/")

# --------------------------------------------------------------------------
# additional linker flags required by openms_gui
if (MSVC)
  ## treat warning of unused function parameter as error, similar to -Werror=unused-variable on GCC
  target_compile_options(OpenMS_GUI PRIVATE "/we4100")
  ## treat warning of unused local variable as error, similar to -Werror=unused-variable on GCC
  target_compile_options(OpenMS_GUI PRIVATE "/we4189")

  set (GUI_lnk_flags "/FORCE:MULTIPLE /INCREMENTAL:NO /ignore:4006 /ignore:4088")
  set_property(TARGET OpenMS_GUI APPEND PROPERTY LINK_FLAGS ${GUI_lnk_flags}) ## allow multiple definitions of symbols (e.g. from template instanciations or STL-derived classes)
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
  ## add icons to TOPPView, INIFileEditor, etc.
  if (MSVC AND EXISTS ${resource_file})
    message(STATUS "Setting resource file ${resource_file} for ${i}")
    add_executable(${i} ${GUI_DIR}/${i}.cpp ${resource_file})
    install_tool(${i})
  elseif (APPLE AND EXISTS ${resource_dir})
    add_mac_app_bundle(${i})
  else()
    if (MSVC)
      message(STATUS "No resource file (${resource_file}) found for ${i}. No icon will be embedded.")
    elseif (APPLE)
      message(STATUS "No icon file (${i}.icns) found in ${resource_dir}. Will not build ${i} as app bundle.")
    endif()
    add_executable(${i} ${GUI_DIR}/${i}.cpp)
    install_tool(${i})
  endif()

  ## append visual lib as dependency for GUI tools
  target_link_libraries(${i} OpenMS_GUI)
endforeach(i)

  ## install platform plugins. On Mac this is handled for every GUI app inside its app bundle.
  if(WIN32)
      install(IMPORTED_RUNTIME_ARTIFACTS "Qt5::QWindowsIntegrationPlugin"
              DESTINATION "${INSTALL_PLUGIN_DIR}/platforms"
              RUNTIME_DEPENDENCY_SET OPENMS_DEPS
              COMPONENT Dependencies)
      file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
        "[Paths]\nPlugins = ${QT_PLUGIN_PATH_REL_TO_BIN}\n")
      install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
        DESTINATION ${INSTALL_BIN_DIR}
        COMPONENT QtPlatformPlugin)
  endif()

add_custom_target(GUI)
add_dependencies(GUI ${GUI_executables})

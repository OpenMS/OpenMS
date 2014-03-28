# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

# --------------------------------------------------------------------------
# before we do anything we try to find Qt since we rely heavily on its
# availability
find_package(Qt4 REQUIRED QtGui QtOpenGL QtSVG QtWebKit QtNetwork ${phonon_extra})
if (NOT QT4_FOUND)
  message(STATUS "Qt4 modules not found!")
	message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
endif()
include(${QT_USE_FILE})

# --------------------------------------------------------------------------
# include custom macros for openms_gui lib
include(qt_wrap_ui.cmake)
include(add_mac_bundle.cmake)

#------------------------------------------------------------------------------
# we need to handle boost > 1.52 when moc-ing so we explicitely find boost here
find_boost()

if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost was not found!")
endif()

# --------------------------------------------------------------------------
# big include file for headers and cpp files
include (${PROJECT_SOURCE_DIR}/includes.cmake)

# --------------------------------------------------------------------------
# Expose OpenMSGUI's public includes (including transitive dependencies)
set(OpenMS_GUI_INTERNAL_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include
                                            ${PROJECT_BINARY_DIR}/include)

set(OpenMS_GUI_EXTERNAL_INCLUDE_DIRECTORIES ${OpenMS_INCLUDE_DIRECTORIES})

set(OpenMS_GUI_INCLUDE_DIRECTORIES ${OpenMS_GUI_INTERNAL_INCLUDE_DIRECTORIES}
                                   ${OpenMS_GUI_EXTERNAL_INCLUDE_DIRECTORIES}
    CACHE INTERNAL "OpenMS_GUI include directories" FORCE)

# add OpenMS include directories
include_directories(${OpenMS_GUI_INTERNAL_INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${OpenMS_GUI_EXTERNAL_INCLUDE_DIRECTORIES})

# --------------------------------------------------------------------------
# OpenMS_GUI lib
oms_add_library(OpenMS_GUI OpenMSVisual_sources)

# we also want to install OpenMS_GUI
install_library(OpenMS_GUI)
install_headers(${OpenMSVisual_sources_h} OpenMS_GUI)

#------------------------------------------------------------------------------
# copy dll to test/doc bin folder on MSVC systems
copy_dll_to_extern_bin(OpenMS_GUI)

# --------------------------------------------------------------------------
# OpenMS_GUI's additional dependencies
# --------------------------------------------------------------------------

# obsolete when CMake MinRequiredVersion becomes >= 2.8.5
if ("${CMAKE_VERSION}" VERSION_LESS "2.8.5" AND NOT(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")) ##
	set(phonon_extra QtPhonon)
	message(STATUS "Trying to find Phonon explicitly... If you run into trouble with Qt Phonon, try to disable it in build system")
else()
	set(phonon_extra)
	message(STATUS "Not trying to find Phonon explicitly... If you run into trouble with Qt Phonon, try to enable it in build system")
endif()

# --------------------------------------------------------------------------
# additional linker flags required by openms_gui
if (MSVC)
	set (GUI_lnk_flags "/FORCE:MULTIPLE /INCREMENTAL:NO /ignore:4006 /ignore:4088")
	set_target_properties(OpenMS_GUI PROPERTIES LINK_FLAGS ${GUI_lnk_flags}) ## allow multiple definitions of symbols (e.g. from template instanciations or STL-derived classes)
endif()

target_link_libraries(OpenMS_GUI OpenMS ${QT_LIBRARIES})

# --------------------------------------------------------------------------
# link libraries for GUI Tools
set(OPENMS_GUI_LIBRARIES OpenMS_GUI OpenMS ${QT_LIBRARIES} CACHE INTERNAL "OpenMS GUI libraries" FORCE)

# --------------------------------------------------------------------------
# GUI tools
set(GUI_DIR source/VISUAL/APPLICATIONS/GUITOOLS)
include(${GUI_DIR}/executables.cmake)
foreach(i ${GUI_executables})
	set(resource_file ${PROJECT_SOURCE_DIR}/${GUI_DIR}/${i}.rc)
	set(resource_dir ${PROJECT_SOURCE_DIR}/${GUI_DIR}/${i}-resources/)
	## add icons to TOPPView and INIFileEditor
	if (MSVC AND EXISTS ${resource_file})
		Message(STATUS "Setting resource file ${resource_file} for ${i}")
		add_executable(${i} ${GUI_DIR}/${i}.cpp ${resource_file})
	elseif(APPLE AND EXISTS ${resource_dir})
		add_mac_app_bundle(${i})
	else()
    if(APPLE)
      message(STATUS "No icon file (${i}.icns) found in ${resource_dir}. Will not build ${i} as app bundle.")
    else()
  		message(STATUS "No resource file (${resource_file}) found for ${i}. No icon will be embedded.")
    endif()
		add_executable(${i} ${GUI_DIR}/${i}.cpp)
	endif()

	## append visual lib as dependency for GUI tools
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
	if (OPENMP_FOUND AND NOT MSVC)
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()

  # we also want to install each GUI tool
  install_tool(${i})
endforeach(i)
add_custom_target(GUI)
add_dependencies(GUI ${GUI_executables})

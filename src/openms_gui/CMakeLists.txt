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
#include(qt_wrap_ui.cmake)
include(qt_wrap_ui.cmake)

# --------------------------------------------------------------------------
# big include file for headers and C files
include (${PROJECT_SOURCE_DIR}/includes.cmake)

# --------------------------------------------------------------------------
# Expose OpenMSGUI's public includes (including transitive dependencies)
set(OPENMS_GUI_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include
                                   ${PROJECT_BINARY_DIR}/include
                                   ${OPENMS_INCLUDE_DIRECTORIES}
    CACHE INTERNAL "OpenMS_GUI include directories" FORCE)

# add OpenMS include directories
include_directories(${OPENMS_GUI_INCLUDE_DIRECTORIES})

# --------------------------------------------------------------------------
# OpenMS_GUI lib
add_library(OpenMS_GUI ${OpenMSVisual_sources})

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

## link libraries for GUI Tools
set(OPENMS_GUI_LIBRARIES optimized OpenMS_GUI${CMAKE_RELEASE_POSTFIX} debug OpenMS_GUI${CMAKE_DEBUG_POSTFIX} OpenMS ${QT_LIBRARIES}
    CACHE INTERNAL "OpenMS GUI libraries" FORCE)
    
##### GUI tools ########
set(GUI_DIR source/VISUAL/APPLICATIONS/GUITOOLS)
include(${GUI_DIR}/executables.cmake)
foreach(i ${GUI_executables})
	set(resource_file ${PROJECT_SOURCE_DIR}/${GUI_DIR}/${i}.rc)
	set(resource_dir ${PROJECT_SOURCE_DIR}/${GUI_DIR}/${i}-resources/)
	## add icons to TOPPView and INIFileEditor
	if (MSVC AND EXISTS ${resource_file})
		Message(STATUS "Setting resource file ${resource_file} for ${i}")
		add_executable(${i} ${GUI_DIR}/${i}.C ${resource_file})
	elseif(APPLE AND EXISTS ${resource_dir})
		# the icon file
		set(ICON_FILE_PATH "${PROJECT_SOURCE_DIR}/source/VISUAL/APPLICATIONS/GUITOOLS/${i}-resources/${i}.icns")
    set(INFO_PLIST_TEMPLATE "${PROJECT_SOURCE_DIR}/source/VISUAL/APPLICATIONS/GUITOOLS/${i}-resources/${i}.plist.in")
		get_filename_component(ICON_FILE_NAME "${ICON_FILE_PATH}" NAME)

		# we also need the icns in the app
		add_executable(
			${i}
			MACOSX_BUNDLE
			${GUI_DIR}/${i}.C
			${ICON_FILE_PATH})

		set_target_properties(${i} PROPERTIES
      # we want our own info.plist template
      MACOSX_BUNDLE_INFO_PLIST "${INFO_PLIST_TEMPLATE}"
			MACOSX_BUNDLE_INFO_STRING "${PROJECT_NAME} Version ${CF_OPENMS_PACKAGE_VERSION}, Copyright 2013 The OpenMS Team."
			MACOSX_BUNDLE_ICON_FILE ${ICON_FILE_NAME}
			MACOSX_BUNDLE_GUI_IDENTIFIER "de.openms.${i}"
			MACOSX_BUNDLE_LONG_VERSION_STRING "${PROJECT_NAME} Version ${CF_OPENMS_PACKAGE_VERSION}"
			MACOSX_BUNDLE_BUNDLE_NAME ${i}
			MACOSX_BUNDLE_SHORT_VERSION_STRING ${CF_OPENMS_PACKAGE_VERSION}
			MACOSX_BUNDLE_BUNDLE_VERSION ${CF_OPENMS_PACKAGE_VERSION}
			MACOSX_BUNDLE_COPYRIGHT "Copyright 2013, The OpenMS Team. All Rights Reserved."
		)

		set_source_files_properties(${ICON_FILE_PATH} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
	else()
    if(APPLE)
      message(STATUS "No icon file (${i}.icns) found in ${resource_dir}. Will not build ${i} as app bundle.")
    else()
  		message(STATUS "No resource file (${resource_file}) found for ${i}. No icon will be embedded.")
    endif()
		add_executable(${i} ${GUI_DIR}/${i}.C)
	endif()

	## append visual lib as dependency for GUI tools
	target_link_libraries(${i} ${OPENMS_GUI_LIBRARIES})
	if (OPENMP_FOUND AND NOT MSVC)
		set_target_properties(${i} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	endif()
endforeach(i)
add_custom_target(GUI)
add_dependencies(GUI ${GUI_executables})
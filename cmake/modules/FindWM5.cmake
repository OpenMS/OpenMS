# - Geometric Tools - WildMagic5 find module
#
# The follwoing variables are optionally searched for defaults
#  WM5_DIR:            Base directory of WM5 installation tree.
#
# The following are set after configuration is done:
#  WM5_FOUND
#  WM5_INCLUDE_DIRS
#  WM5_LIBRARIES
#  WM5_LINK_DIRECTORIES
#  WM5_DEFINITIONS
#

#=============================================================================
# Copyright 2013 The OpenMS Team -- Eberhard Karls University Tuebingen,
#                ETH Zurich, and Freie Universitaet Berlin.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file LICENSE for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

#------------------------------------------------------------------------------
# List of WM5 components relevant for OpenMS
set(WM5_LIBS Wm5Core Wm5Mathematics)

#------------------------------------------------------------------------------
# find WM5 include path
find_path(_WM5_INCLUDE_DIR
          NAMES
            Wm5Core.h
          HINTS
            ${WM5_DIR}/include
          PATH_SUFFIXES WildMagic
)

#------------------------------------------------------------------------------
# find the relevant WM5 libraries
set(_WM5_LIBRARIES "")
set(_ALL_WM5_LIBS_FOUND TRUE)

# find all components
foreach(_lib ${WM5_LIBS})
	string(TOUPPER WM5_${_lib} _WM5LIB)
  message(STATUS "Searching for lib: ${_lib}")

  find_library(${_WM5LIB}_LIBRARY_RELEASE
    ${_lib}
    HINTS ${WM5_DIR}/lib
  )

  find_library(${_WM5LIB}_LIBRARY_DEBUG
    ${_lib}d
    HINTS ${WM5_DIR}/lib
  )

  # select debug or release depending on current config
  include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)
  select_library_configurations(${_WM5LIB})

  if(${_WM5LIB}_LIBRARY)
    message(STATUS "Found ${_lib}: ${${_WM5LIB}_LIBRARY}")
  else(${_WM5LIB}_LIBRARY)
    message(STATUS "WARNING: ${_lib} was not found!")
    set(_ALL_WM5_LIBS_FOUND FALSE)
  endif(${_WM5LIB}_LIBRARY)

  set(_WM5_LIBRARIES ${WM5_LIBRARIES} ${${_WM5LIB}_LIBRARY})
endforeach(_lib ${WM5_LIBS})

#------------------------------------------------------------------------------
# mark as found if we have a found the include path and all libraries
if(_WM5_INCLUDE_DIR AND _ALL_WM5_LIBS_FOUND)
  # we assume we found the correct version
  set(WM5_PROPER_VERSION_FOUND TRUE)
  set(WM5_FOUND TRUE)
else()
  set(WM5_FOUND FALSE)
endif()

#------------------------------------------------------------------------------
# check if we have a static Wm5 lib
if(MSVC)
	# we try the optimized version of WM5_CORE to win32 path
	list(GET WM5_WM5CORE_LIBRARY 1 _wm5_opt_core)
	file(TO_NATIVE_PATH "${_wm5_opt_core}" _native_wm5_opt_core)

	# find dumpbin and translate to win32 path
	find_program( _dumpbin_executable "dumpbin.exe" PATHS ENV PATH )
	file(TO_NATIVE_PATH "${_dumpbin_executable}" _native_dumpbin_executable)

	if(NOT _dumpbin_executable)
		message(STATUS      "Failed to find dumpbin.exe.")
		message(FATAL_ERROR "Please make sure dumpbin is in your PATH and that you use a Visual Studio Command Line.")
	endif()

	# execute dumpbin /exports -> shared libs will contain "Exports"
	execute_process(COMMAND ${_native_dumpbin_executable} "/exports" ${_native_wm5_opt_core}
		OUTPUT_VARIABLE _wm5_dumpbin_stdout
		ERROR_VARIABLE _wm5_dumpbin_stderr
		RESULT_VARIABLE  _wm5_dumpbin_result
	)

	# dumpbin call failed -> print diagnostics for debugging
	if(NOT _wm5_dumpbin_result EQUAL 0)
                message(STATUS "Failed to execute ${_native_dumpbin_executable} to assess Wm5 library type.")
                message(STATUS "dumpbin returned: ${_wm5_dumpbin_result}")
                message(STATUS "dumpbin stdout: ${_wm5_dumpbin_stdout}")
                message(FATAL_ERROR "dumpbin stderr: ${_wm5_dumpbin_stderr}")
                message(STATUS $ENV{PATH})
	endif()

	# check if result contains "Exports"
	string(FIND "Exports" ${_wm5_dumpbin_stdout} _export_match)
	if(_export_match EQUAL -1)
		set(WM5_IS_STATIC TRUE)
	else()
		set(WM5_IS_STATIC FALSE)
	endif()
else()
	get_filename_component(wm5core_lib_ext "${WM5_WM5CORE_LIBRARY}" EXT)
	message(STATUS "wm5core_lib_ext: ${wm5core_lib_ext}")
	if("${wm5core_lib_ext}" STREQUAL "${CMAKE_STATIC_LIBRARY_SUFFIX}")
		set(WM5_IS_STATIC TRUE)
	endif()
endif()

#------------------------------------------------------------------------------
# print warning/error if NOT FIND_QUIETLY and NOT found
if(NOT WM5_FOUND AND NOT WM5_FIND_QUIETLY)
  if(WM5_FIND_REQUIRED)
    message(FATAL_ERROR "WM5 required but wasn't found completely.")
  else(WM5_FIND_REQUIRED)
    message(STATUS "WARNING: WM5 was not found.")
  endif(WM5_FIND_REQUIRED)
endif(NOT WM5_FOUND AND NOT WM5_FIND_QUIETLY)

#------------------------------------------------------------------------------
# add definitions for dll linking on windows
set(_WM5_DEFINITIONS "")
if(MSVC AND NOT WM5_IS_STATIC)
  set(_WM5_DEFINITIONS "/D WM5_CORE_DLL_IMPORT /D WM5_MATHEMATICS_DLL_IMPORT")
endif()

#------------------------------------------------------------------------------
# handle standard args
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(WM5
                                  DEFAULT_MSG
                                  _WM5_INCLUDE_DIR
                                  _WM5_LIBRARIES
                                  WM5_PROPER_VERSION_FOUND)

#------------------------------------------------------------------------------
# export found libraries
if(WM5_FOUND)
  set(WM5_INCLUDE_DIRS ${_WM5_INCLUDE_DIR})
  set(WM5_LIBRARIES ${_WM5_LIBRARIES})
  set(WM5_DEFINITIONS ${_WM5_DEFINITIONS})
endif(WM5_FOUND)

#------------------------------------------------------------------------------
# hide internal entries
mark_as_advanced(_WM5_LIBRARIES
                 _WM5_INCLUDE_DIR
                 _WM5_DEFINITIONS)

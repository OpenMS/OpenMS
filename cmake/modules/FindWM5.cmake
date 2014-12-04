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
# find WM5 include path
find_path(WM5_INCLUDE_DIR
          NAMES
            Wm5Core.h
          HINTS
            ${WM5_DIR}/include
          PATH_SUFFIXES WildMagic libwildmagic
)

#------------------------------------------------------------------------------
# find the relevant WM5 libraries
set(WM5_LIBRARIES "")
set(_ALL_WM5_LIBS_FOUND TRUE)

include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)

# helper macro to find specific coind sub-libraries
macro(_wm5_find_lib _libname)
  if(NOT WM5_${_libname})
    # find release version
    find_library(WM5_${_libname}_LIBRARY_RELEASE
      NAMES ${_libname}
      HINTS ${WM5_DIR}/lib
    )
    # .. and debug version
    find_library(WM5_${_libname}_LIBRARY_DEBUG
      NAMES ${_libname}d
      HINTS ${WM5_DIR}/lib
    )

    # create final library to be exported
    select_library_configurations(WM5_${_libname})
  endif()
endmacro()

#------------------------------------------------------------------------------
# find Wm5 libs relevant for OpenMS
_wm5_find_lib("Wm5Core")
_wm5_find_lib("Wm5Mathematics")

#------------------------------------------------------------------------------
# handle standard args
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(WM5 DEFAULT_MSG
  WM5_INCLUDE_DIR
  WM5_Wm5Core_LIBRARY
  WM5_Wm5Mathematics_LIBRARY)

# export the libraries and handle dynamic linkage stuff
if(WM5_FOUND)
  set(WM5_INCLUDE_DIRS ${WM5_INCLUDE_DIR})
  set(WM5_LIBRARIES
    ${WM5_Wm5Core_LIBRARY}
    ${WM5_Wm5Mathematics_LIBRARY})

  #------------------------------------------------------------------------------
  # check if we have a static Wm5 lib
  if(MSVC AND NOT DEFINED WM5_IS_STATIC)
    # we try the optimized version of WM5_CORE to win32 path
    list(GET WM5_Wm5Core_LIBRARY 1 _wm5_opt_core)
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
    get_filename_component(wm5core_lib_ext "${WM5_Wm5Core_LIBRARY}" EXT)
    if("${wm5core_lib_ext}" STREQUAL "${CMAKE_STATIC_LIBRARY_SUFFIX}")
      set(WM5_IS_STATIC TRUE)
    endif()
  endif()

  #------------------------------------------------------------------------------
  # add definitions for dll linking on windows
  set(WM5_DEFINITIONS "")
  if(MSVC AND NOT WM5_IS_STATIC)
    set(WM5_DEFINITIONS "/D WM5_CORE_DLL_IMPORT /D WM5_MATHEMATICS_DLL_IMPORT")
  endif()
endif()

#------------------------------------------------------------------------------
# hide internal entries
mark_as_advanced(_WM5_LIBRARIES
                 _WM5_INCLUDE_DIR
                 _WM5_DEFINITIONS
                 WM5_Wm5Core_LIBRARIES
                 WM5_Wm5Mathematics_LIBRARIES
                 WM5_IS_STATIC)

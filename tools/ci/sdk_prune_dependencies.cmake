# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
#
# Reduce the bundled runtime dependencies of the developer SDK to what the SDK's
# own libraries load. Called by tools/ci/package_sdk.sh:
#
#   cmake -DLIB_DIR=<SDK lib dir> "-DROOTS=<lib>;<lib>;..." -P sdk_prune_dependencies.cmake
#
# The install component Dependencies is the runtime dependency set of every
# installed target (cmake/package_general.cmake), so it also holds what only the
# GUI library and the applications need: Qt and its plugins on Windows and macOS,
# among others. The SDK ships the core and CLI layers only, so this follows the
# dependencies of ROOTS (those layers' libraries) the way the loader does, and
# removes every shared library, framework or plugin below LIB_DIR outside that
# closure. Static libraries, import libraries and the CMake package files are
# left alone.
#
# A library that exists in LIB_DIR but cannot be resolved stops the script
# instead of being removed: the closure is only trusted when it is complete.
cmake_minimum_required(VERSION 3.24)

foreach(_var IN ITEMS LIB_DIR ROOTS)
  if(NOT DEFINED ${_var} OR "${${_var}}" STREQUAL "")
    message(FATAL_ERROR "sdk_prune_dependencies.cmake: ${_var} is required")
  endif()
endforeach()
file(REAL_PATH "${LIB_DIR}" LIB_DIR)

# In script mode, CMake has to be told how to read the binaries of this platform
# (an install script gets this from the configured project).
if(CMAKE_HOST_WIN32)
  set(CMAKE_GET_RUNTIME_DEPENDENCIES_PLATFORM "windows+pe")
  set(CMAKE_GET_RUNTIME_DEPENDENCIES_TOOL "dumpbin")
elseif(CMAKE_HOST_APPLE)
  set(CMAKE_GET_RUNTIME_DEPENDENCIES_PLATFORM "macos+macho")
  set(CMAKE_GET_RUNTIME_DEPENDENCIES_TOOL "otool")
else()
  set(CMAKE_GET_RUNTIME_DEPENDENCIES_PLATFORM "linux+elf")
  set(CMAKE_GET_RUNTIME_DEPENDENCIES_TOOL "objdump")
endif()

# system libraries are neither bundled nor followed (cf. cmake/package_general.cmake)
set(_pre_exclude "^api-ms-" "^ext-ms-" "^hvsi" "^pdmutilities")
set(_post_exclude
  ".*[\\/][Ss][Yy][Ss][Tt][Ee][Mm]32[\\/].*"
  "^/usr/lib/" "^/System/" "^/lib/" "^/lib64/")

# A dependency that different libraries resolve to different files is reported as
# a conflict instead of a resolved dependency, and its own dependencies are not
# followed. The SDK's copy is what its libraries load, so such a copy becomes a
# root of the next pass, until no conflict involves the SDK.
set(_roots ${ROOTS})
set(_pass 0)
while(TRUE)
  math(EXPR _pass "${_pass} + 1")
  if(_pass GREATER 10)
    message(FATAL_ERROR "Conflicting dependencies keep appearing after ${_pass} passes: ${_conflict_FILENAMES}")
  endif()
  file(GET_RUNTIME_DEPENDENCIES
    LIBRARIES ${_roots}
    RESOLVED_DEPENDENCIES_VAR _resolved
    UNRESOLVED_DEPENDENCIES_VAR _unresolved
    CONFLICTING_DEPENDENCIES_PREFIX _conflict
    DIRECTORIES "${LIB_DIR}"
    PRE_EXCLUDE_REGEXES ${_pre_exclude}
    POST_EXCLUDE_REGEXES ${_post_exclude})
  set(_new_roots)
  foreach(_name IN LISTS _conflict_FILENAMES)
    message(STATUS "found in several places: ${_name}: ${_conflict_${_name}}")
    foreach(_path IN LISTS _conflict_${_name})
      file(REAL_PATH "${_path}" _real)
      string(FIND "${_real}" "${LIB_DIR}/" _pos)
      if(_pos EQUAL 0 AND NOT _real IN_LIST _roots)
        list(APPEND _new_roots "${_real}")
      endif()
    endforeach()
  endforeach()
  if(NOT _new_roots)
    break()
  endif()
  list(APPEND _roots ${_new_roots})
endwhile()

# an unresolved name that the SDK does contain means the closure is incomplete
set(_missing)
foreach(_name IN LISTS _unresolved)
  get_filename_component(_base "${_name}" NAME)
  file(GLOB_RECURSE _present "${LIB_DIR}/${_base}")
  if(_present)
    list(APPEND _missing "${_name}")
  else()
    message(STATUS "unresolved, not part of the SDK: ${_name}")
  endif()
endforeach()
if(_missing)
  message(FATAL_ERROR "Could not resolve these libraries although the SDK contains them, "
                      "so nothing was removed: ${_missing}")
endif()

set(_keep)
foreach(_file IN LISTS _roots _resolved)
  file(REAL_PATH "${_file}" _real)
  list(APPEND _keep "${_real}")
endforeach()

# frameworks as a whole: kept when the closure uses their binary
file(GLOB_RECURSE _entries LIST_DIRECTORIES true "${LIB_DIR}/*")
foreach(_entry IN LISTS _entries)
  if(IS_DIRECTORY "${_entry}" AND _entry MATCHES "\\.framework$"
     AND NOT _entry MATCHES "\\.framework/.*\\.framework$")
    set(_used FALSE)
    foreach(_kept IN LISTS _keep)
      string(FIND "${_kept}" "${_entry}/" _pos)
      if(_pos EQUAL 0)
        set(_used TRUE)
        break()
      endif()
    endforeach()
    if(NOT _used)
      message(STATUS "removing ${_entry}")
      file(REMOVE_RECURSE "${_entry}")
    endif()
  endif()
endforeach()

# shared libraries and the symbolic links naming them: decided before anything is
# removed, so a link is judged by the file it points to, not by whether that file
# has already gone
file(GLOB_RECURSE _files LIST_DIRECTORIES false "${LIB_DIR}/*")
set(_remove)
foreach(_file IN LISTS _files)
  if(_file MATCHES "\\.framework/" OR NOT EXISTS "${_file}")
    continue()
  endif()
  if(NOT _file MATCHES "\\.(so(\\.[0-9]+)*|dylib|dll)$")
    continue()
  endif()
  file(REAL_PATH "${_file}" _real)
  if(NOT _real IN_LIST _keep)
    list(APPEND _remove "${_file}")
  endif()
endforeach()
foreach(_file IN LISTS _remove)
  message(STATUS "removing ${_file}")
  file(REMOVE "${_file}")
endforeach()

# directories left empty (e.g. the Qt plugin directory), deepest first
file(GLOB_RECURSE _dirs LIST_DIRECTORIES true "${LIB_DIR}/*")
list(SORT _dirs ORDER DESCENDING)
foreach(_dir IN LISTS _dirs)
  if(IS_DIRECTORY "${_dir}" AND NOT IS_SYMLINK "${_dir}")
    file(GLOB _content "${_dir}/*")
    if(NOT _content)
      file(REMOVE_RECURSE "${_dir}")
    endif()
  endif()
endforeach()

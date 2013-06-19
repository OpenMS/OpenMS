#  -- Try to find LibSVM
#
#  LibSVM is a Library for Support Vector Machines
#  available at http://www.csie.ntu.edu.tw/~cjlin/libsvm/
#
#  ------------------------------------------------------------------
#
#  -- Library usage example :
#
#  find_package (LibSVM 2.9.0)
#  if (LIBSVM_FOUND)
#     include_directories (${LIBSVM_INCLUDE_DIRS})
#     add_executable (foo foo.cpp)
#     target_link_libraries (foo ${LIBSVM_LIBRARIES})
#  endif ()
#
#  -- Show some debug information :
#
#  set (LIBSVM_DEBUG TRUE)
#  find_package (LibSVM)
#
#  -------------------------------------------------------------------
#
#  -- This module defines :
#
#  LIBSVM_FOUND - the system has LibSVM
#  LIBSVM_INCLUDE_DIR - where to find svm.h
#  LIBSVM_INCLUDE_DIRS libsvm includes
#  LIBSVM_LIBRARY - where to find the LibSVM library
#  LIBSVM_LIBRARIES - aditional libraries
#  LIBSVM_VERSION - version
#  LIBSVM_MAJOR_VERSION - major version
#  LIBSVM_MINOR_VERSION - minor version
#  LIBSVM_SUBMINOR_VERSION - subminor version
#  LIBSVM_ROOT_DIR - root dir (ex. /usr/local)
#
#  ------------------------------------------------------------------
#
#  Copyright (c) 2010 Julien Schueller <schueller at phimeca dot com>
#
#  Redistribution and use is allowed according to the terms of the New
#  BSD license.
#  For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

# check for math header
if (NOT HAVE_MATH_H)
  include (CheckIncludeFile)
  check_include_file (math.h HAVE_MATH_H)
endif ()


# set LIBSVM_INCLUDE_DIR
if (NOT LIBSVM_INCLUDE_DIR)
  find_path (LIBSVM_INCLUDE_DIR
    NAMES
      svm.h
    PATHS
      ${LIBSVM_ROOT_DIR}/include
    PATH_SUFFIXES
      libsvm
      libsvm-3.1/libsvm
    DOC
      "LibSVM include directory"
  )
endif ()


# set LIBSVM_INCLUDE_DIR
if (NOT LIBSVM_INCLUDE_DIRS)
  set (LIBSVM_INCLUDE_DIRS ${LIBSVM_INCLUDE_DIR})
endif ()


# version
if (NOT LIBSVM_VERSION)
  set (LIBSVM_MAJOR_VERSION 0)
  set (LIBSVM_MINOR_VERSION 0)
  set (LIBSVM_SUBMINOR_VERSION 0)
  if (LIBSVM_INCLUDE_DIR)
    # LIBSVM_VERSION macro defined in svm.h since version 2.8.9
    file (STRINGS "${LIBSVM_INCLUDE_DIR}/svm.h" _VERSION_STRING REGEX ".*LIBSVM_VERSION.*")
    if (_VERSION_STRING)
      string (REGEX REPLACE ".*_VERSION[ ]+([0-9]+)" "\\1" _VERSION_NUMBER "${_VERSION_STRING}")
      math (EXPR LIBSVM_MAJOR_VERSION "${_VERSION_NUMBER} / 100")
      math (EXPR LIBSVM_MINOR_VERSION "(${_VERSION_NUMBER} % 100 ) / 10")
      math (EXPR LIBSVM_SUBMINOR_VERSION "${_VERSION_NUMBER} % 10")
    endif ()
  endif ()
  set (LIBSVM_VERSION "${LIBSVM_MAJOR_VERSION}.${LIBSVM_MINOR_VERSION}.${LIBSVM_SUBMINOR_VERSION}")
endif ()



# check version
set (_LIBSVM_VERSION_MATCH TRUE)
set (_REQUIRED_VERSION "${LibSVM_FIND_VERSION_MAJOR}.${LibSVM_FIND_VERSION_MINOR}.${LibSVM_FIND_VERSION_PATCH}")
if (LibSVM_FIND_VERSION AND _VERSION_STRING)
  if (LibSVM_FIND_VERSION_EXACT)
    if ("${_REQUIRED_VERSION}" VERSION_EQUAL "${LIBSVM_VERSION}")
    else()
      set (_LIBSVM_VERSION_MATCH FALSE)
    endif ()
  else ()
    if ("${_REQUIRED_VERSION}" VERSION_GREATER "${LIBSVM_VERSION}")
      set (_LIBSVM_VERSION_MATCH FALSE)
    endif ()
  endif ()
endif ()


# set LIBSVM_LIBRARY
if(NOT LIBSVM_LIBRARY)
  find_library (LIBSVM_LIBRARY
    NAMES
      svm
    PATHS
      ${LIBSVM_ROOT_DIR}/lib
    DOC
      "LibSVM library location"
  )
endif ()


# set LIBSVM_LIBRARIES
if (NOT LIBSVM_LIBRARIES)
  set (LIBSVM_LIBRARIES ${LIBSVM_LIBRARY})

  # link with math library on unix
  if (UNIX)
    list (APPEND LIBSVM_LIBRARIES "-lm")
  endif ()
endif ()


# root dir
if (NOT LIBSVM_ROOT_DIR)
  # try to guess root dir from include dir
  if (LIBSVM_INCLUDE_DIR)
    string (REGEX REPLACE "(.*)/include.*" "\\1" LIBSVM_ROOT_DIR ${LIBSVM_INCLUDE_DIR})

  # try to guess root dir from library dir
  elseif (LIBSVM_LIBRARY)
    string (REGEX REPLACE "(.*)/lib[/|32|64].*" "\\1" LIBSVM_ROOT_DIR ${LIBSVM_LIBRARY})
  endif ()
endif ()


# debug
if (LIBSVM_DEBUG)
  message (STATUS "LIBSVM_LIBRARY: ${LIBSVM_LIBRARY}")
  message (STATUS "LIBSVM_LIBRARIES: ${LIBSVM_LIBRARIES}")
  message (STATUS "LIBSVM_INCLUDE_DIR: ${LIBSVM_INCLUDE_DIR}")
  message (STATUS "LIBSVM_INCLUDE_DIRS: ${LIBSVM_INCLUDE_DIRS}")
  message (STATUS "LIBSVM_ROOT_DIR: ${LIBSVM_ROOT_DIR}")
  message (STATUS "LIBSVM_VERSION: ${LIBSVM_VERSION}")
  message (STATUS "LIBSVM_MAJOR_VERSION: ${LIBSVM_MAJOR_VERSION}")
  message (STATUS "LIBSVM_MINOR_VERSION: ${LIBSVM_MINOR_VERSION}")
  message (STATUS "LIBSVM_SUBMINOR_VERSION: ${LIBSVM_SUBMINOR_VERSION}")
endif ()


# handle REQUIRED and QUIET options
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (LibSVM DEFAULT_MSG LIBSVM_LIBRARY
  _LIBSVM_VERSION_MATCH
  HAVE_MATH_H
  LIBSVM_INCLUDE_DIR
  LIBSVM_INCLUDE_DIRS
  LIBSVM_LIBRARIES
  LIBSVM_ROOT_DIR
  LIBSVM_VERSION
)


mark_as_advanced (
  LIBSVM_LIBRARY
  LIBSVM_LIBRARIES
  LIBSVM_INCLUDE_DIR
  LIBSVM_INCLUDE_DIRS
  LIBSVM_ROOT_DIR
  LIBSVM_VERSION
)

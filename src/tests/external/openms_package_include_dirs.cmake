# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause

### Helper for the installed-consumer tests of src/tests/external. A plain
### consumer project needs none of this: it links the imported targets and CMake
### resolves their include directories itself.
###
### openms_package_include_dirs(<target> <out_var>)
###
### The include directories an installed OpenMS package hands a consumer, as
### plain paths.
###
### Reading INTERFACE_INCLUDE_DIRECTORIES is not enough. The public headers are
### installed as a CMake header file set (cmake/install_macros.cmake), so the
### generated OpenMSTargets.cmake re-creates that file set on the imported
### target with
###
###   target_sources(OpenMS::OpenMS INTERFACE FILE_SET "HEADERS" TYPE "HEADERS"
###                  BASE_DIRS "${_IMPORT_PREFIX}/include" FILES ...)
###
### and CMake records a file set's base directories in
### INTERFACE_INCLUDE_DIRECTORIES wrapped in $<BUILD_INTERFACE:...>. That
### wrapper evaluates to the directory in every context an imported target is
### used in -- it is dropped only where a project exports targets of its own --
### so consumers compile against the right headers either way. Only a
### configure-time reader of the property sees it, and a check that treats the
### raw entry as a path then quietly stops checking anything.
function(openms_package_include_dirs target out_var)
  set(_dirs)

  get_target_property(_includes ${target} INTERFACE_INCLUDE_DIRECTORIES)
  if(_includes)
    foreach(_include IN LISTS _includes)
      ## ".*" rather than ".+": openms_add_library() writes a
      ## "$<BUILD_INTERFACE:${EXTERNAL_INCLUDES}>" that is empty for a library
      ## without external includes, which names no directory at all.
      string(REGEX REPLACE "^\\$<BUILD_INTERFACE:(.*)>$" "\\1" _include "${_include}")
      if(_include MATCHES "\\$<")
        message(FATAL_ERROR
          "an include directory of ${target} is the generator expression '${_include}', "
          "which openms_package_include_dirs() cannot reduce to a path. Teach it about "
          "that expression rather than leaving the checks built on it passing on a "
          "string that is not a directory.")
      endif()
      if(NOT _include STREQUAL "")
        list(APPEND _dirs "${_include}")
      endif()
    endforeach()
  endif()

  ## The same directories as the file sets themselves report them, so that these
  ## checks keep seeing the package's headers no matter how CMake chooses to
  ## surface a file set's base directories in the include directories.
  get_target_property(_header_sets ${target} INTERFACE_HEADER_SETS)
  if(_header_sets)
    foreach(_header_set IN LISTS _header_sets)
      ## The default file set (name and type HEADERS) uses HEADER_DIRS;
      ## HEADER_DIRS_<NAME> is an error for it.
      if(_header_set STREQUAL "HEADERS")
        get_target_property(_base_dirs ${target} HEADER_DIRS)
      else()
        get_target_property(_base_dirs ${target} HEADER_DIRS_${_header_set})
      endif()
      if(_base_dirs)
        list(APPEND _dirs ${_base_dirs})
      endif()
    endforeach()
  endif()

  if(_dirs)
    list(REMOVE_DUPLICATES _dirs)
  endif()
  set(${out_var} "${_dirs}" PARENT_SCOPE)
endfunction()

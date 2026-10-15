# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

# Removes empty entries from the RUNPATH (or RPATH) of every ELF file in the CPack
# staging directory. CPack runs this script (CPACK_PRE_BUILD_SCRIPTS, set in
# cmake/package_deb.cmake) after it has installed the files and before it builds
# the package; it is not included at configure time.
#
# The dynamic loader reads an empty entry as the current working directory, so no
# packaged file may carry one. They appear when the binaries were linked under a
# different configuration than the one being packaged, which is what the release
# workflow does: it builds and tests with PACKAGE_TYPE=none, then reconfigures with
# the package type and runs cpack without relinking. CMake pads a build-tree RPATH
# with ':' to the length of the install RPATH, and the install RPATH of
# PACKAGE_TYPE=none also names the vcpkg directory of the build tree (top-level
# CMakeLists.txt), so it is far longer than that of a package. When installing,
# CMake replaces only the padding it now expects and keeps the rest: a library that
# links nothing outside the system, such as libOpenSwathAlgo, was packaged with
# $ORIGIN/../lib/ followed by about 90 empty entries. Relinking before packaging
# would avoid the padding, but 'all' also contains the documentation target, which
# is always out of date and would regenerate the whole documentation.

## Without it the glob below would start at the file system root.
if(NOT IS_DIRECTORY "${CPACK_TEMPORARY_DIRECTORY}")
  message(FATAL_ERROR "cpack_clean_runpath.cmake: CPACK_TEMPORARY_DIRECTORY ('${CPACK_TEMPORARY_DIRECTORY}') "
                      "is not a directory; this script is meant to be run by CPack (CPACK_PRE_BUILD_SCRIPTS).")
endif()

cmake_policy(PUSH)
## Unset policies take their OLD behavior in a CPack script: CMP0009 NEW keeps
## GLOB_RECURSE from following symlinks out of the staging directory.
cmake_policy(VERSION 3.24)

## A function, so that none of its variables reach the CPack generator code that
## runs in this scope afterwards.
function(_openms_clean_runpath _staging_dir)
  file(GLOB_RECURSE _staged_files LIST_DIRECTORIES false "${_staging_dir}/*")
  foreach(_file IN LISTS _staged_files)
    ## Editing through a symlink could reach a file outside the staging directory;
    ## a staged target is edited under its own name.
    if(IS_SYMLINK "${_file}")
      continue()
    endif()

    ## READ_ELF only sets the variables it finds something for.
    unset(_not_elf)
    unset(_rpath)
    unset(_runpath)
    file(READ_ELF "${_file}" RPATH _rpath RUNPATH _runpath CAPTURE_ERROR _not_elf)
    if(_not_elf)
      continue()
    endif()
    ## The loader ignores RPATH when RUNPATH is present.
    if(DEFINED _runpath)
      set(_old "${_runpath}")
    elseif(DEFINED _rpath)
      set(_old "${_rpath}")
    else()
      continue()
    endif()

    ## READ_ELF returns the entries as a list; compare them as the loader sees them.
    string(REPLACE ";" ":" _old "${_old}")
    string(REGEX REPLACE "::+" ":" _new "${_old}")
    string(REGEX REPLACE "^:|:$" "" _new "${_new}")
    if(_new STREQUAL _old)
      continue()
    endif()

    if(_new STREQUAL "")
      file(RPATH_REMOVE FILE "${_file}")
    else()
      file(RPATH_SET FILE "${_file}" NEW_RPATH "${_new}")
    endif()
  endforeach()
endfunction()

_openms_clean_runpath("${CPACK_TEMPORARY_DIRECTORY}")

cmake_policy(POP)

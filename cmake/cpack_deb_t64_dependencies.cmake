# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------

# Adds the pre-t64 name of each renamed library package to the Depends of every DEB.
# CPack runs this script (CPACK_POST_BUILD_SCRIPTS, set in cmake/package_deb.cmake)
# after it has written the packages and before it copies them to the build tree; it
# is not included at configure time.
#
# Debian and Ubuntu (from 24.04) renamed the library packages whose ABI involves
# time_t when they moved their 32-bit architectures to a 64-bit time_t: libqt6core6
# became libqt6core6t64, and so on. dpkg-shlibdeps names the packages installed on the
# build host, so a package built on such a release depends on the new names, which
# older releases do not have. On architectures whose time_t was 64-bit already (amd64,
# arm64, ...) the ABI did not change, and the renamed package provides and replaces
# its old name. There this script adds the old name as an alternative, with the same
# version constraint:
#   libqt6core6t64 (>= 6.2.2)  ->  libqt6core6t64 (>= 6.2.2) | libqt6core6 (>= 6.2.2)
# The old names come from the build host's dpkg database, which dpkg-shlibdeps has
# just read. On armhf, where the ABI changed, the renamed packages do not provide
# their old names, and the dependencies stay as they are. A package built on a release
# before the renaming needs nothing: on those architectures the renamed packages
# provide the names it depends on.
#
# Only the control member of the package is replaced, in place; the data member is
# kept byte for byte.

## Without it the working directory below would be created at the file system root.
if(NOT IS_DIRECTORY "${CPACK_PACKAGE_DIRECTORY}")
  message(FATAL_ERROR "cpack_deb_t64_dependencies.cmake: CPACK_PACKAGE_DIRECTORY ('${CPACK_PACKAGE_DIRECTORY}') "
                      "is not a directory; this script is meant to be run by CPack (CPACK_POST_BUILD_SCRIPTS).")
endif()

cmake_policy(PUSH)
cmake_policy(VERSION 3.24)

## The names that the installed package <package> (name:arch) both provides and
## replaces, i.e. the names it was renamed from. Empty if there are none.
function(_openms_deb_renamed_from _dpkg_query _package _out)
  foreach(_field IN ITEMS Provides Replaces)
    execute_process(
      COMMAND "${_dpkg_query}" --show "--showformat=\${${_field}}" "${_package}"
      OUTPUT_VARIABLE _value
      RESULT_VARIABLE _result
      ERROR_VARIABLE _error
    )
    if(NOT _result EQUAL 0)
      message(WARNING "cpack_deb_t64_dependencies.cmake: cannot query ${_package}, so its dependency stays as it is: ${_error}")
      set(${_out} "" PARENT_SCOPE)
      return()
    endif()
    ## Relations are comma-separated; drop the version and architecture qualifiers.
    string(REPLACE "," ";" _entries "${_value}")
    set(_${_field})
    foreach(_entry IN LISTS _entries)
      if(_entry MATCHES "^[ \t\r\n]*([^ \t\r\n(:]+)")
        list(APPEND _${_field} "${CMAKE_MATCH_1}")
      endif()
    endforeach()
  endforeach()

  set(_names)
  foreach(_name IN LISTS _Provides)
    if(_name IN_LIST _Replaces)
      list(APPEND _names "${_name}")
    endif()
  endforeach()
  set(${_out} "${_names}" PARENT_SCOPE)
endfunction()

## A function, so that none of its variables reach the CPack generator code that
## runs in this scope afterwards.
function(_openms_deb_add_t64_alternatives _work_dir)
  find_program(_dpkg_deb dpkg-deb REQUIRED NO_CACHE)
  find_program(_dpkg_query dpkg-query REQUIRED NO_CACHE)
  find_program(_ar ar REQUIRED NO_CACHE)
  find_program(_tar tar REQUIRED NO_CACHE)
  find_program(_gzip gzip REQUIRED NO_CACHE)

  foreach(_package IN LISTS ARGN)
    if(NOT _package MATCHES "\\.deb$")
      continue()
    endif()
    get_filename_component(_package_name "${_package}" NAME)

    file(REMOVE_RECURSE "${_work_dir}")
    file(MAKE_DIRECTORY "${_work_dir}/control")

    execute_process(
      COMMAND "${_dpkg_deb}" --ctrl-tarfile "${_package}"
      OUTPUT_FILE "${_work_dir}/control.tar"
      RESULT_VARIABLE _result
      ERROR_VARIABLE _error
    )
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "Could not read the control archive of ${_package}: ${_error}")
    endif()
    ## The members in their order, so that the rebuilt archive lists them the same way.
    execute_process(
      COMMAND "${_tar}" --list --file "${_work_dir}/control.tar"
      OUTPUT_VARIABLE _control_members
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _result
      ERROR_VARIABLE _error
    )
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "Could not list the control archive of ${_package}: ${_error}")
    endif()
    string(REPLACE "\n" ";" _control_members "${_control_members}")
    execute_process(
      COMMAND "${_tar}" --extract --preserve-permissions --file "${_work_dir}/control.tar" -C "${_work_dir}/control"
      RESULT_VARIABLE _result
      ERROR_VARIABLE _error
    )
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "Could not extract the control archive of ${_package}: ${_error}")
    endif()

    set(_control_file "${_work_dir}/control/control")
    file(READ "${_control_file}" _control)
    if(NOT _control MATCHES "(^|\n)Architecture:[ \t]*([^\r\n]*)")
      message(FATAL_ERROR "No Architecture field in the control file of ${_package}")
    endif()
    set(_arch "${CMAKE_MATCH_2}")
    if(NOT _control MATCHES "(^|\n)Depends:[ \t]*([^\r\n]*)")
      continue()
    endif()
    set(_depends_field "${CMAKE_MATCH_0}")
    set(_depends_prefix "${CMAKE_MATCH_1}")
    set(_depends "${CMAKE_MATCH_2}")

    string(REPLACE "," ";" _dependencies "${_depends}")
    set(_new_dependencies)
    foreach(_dependency IN LISTS _dependencies)
      string(STRIP "${_dependency}" _dependency)
      ## dpkg-shlibdeps writes "name" or "name (op version)"; entries that already have
      ## alternatives are left alone.
      if(_dependency MATCHES "^([a-z0-9][a-z0-9+.-]*t64)([ \t]*\\([^)]*\\))?$")
        set(_version "${CMAKE_MATCH_2}")
        _openms_deb_renamed_from("${_dpkg_query}" "${CMAKE_MATCH_1}:${_arch}" _old_names)
        foreach(_old_name IN LISTS _old_names)
          string(APPEND _dependency " | ${_old_name}${_version}")
        endforeach()
      endif()
      list(APPEND _new_dependencies "${_dependency}")
    endforeach()
    list(JOIN _new_dependencies ", " _new_depends)
    if(_new_depends STREQUAL _depends)
      continue()
    endif()

    string(REPLACE "${_depends_field}" "${_depends_prefix}Depends: ${_new_depends}" _control "${_control}")
    file(WRITE "${_control_file}" "${_control}")

    ## CPack writes the control archive as control.tar.gz. Rebuild it the same way:
    ## same members in the same order, owned by root, dated no later than
    ## SOURCE_DATE_EPOCH if that is set, and compressed without a time stamp.
    set(_mtime)
    if(DEFINED ENV{SOURCE_DATE_EPOCH})
      set(_mtime "--mtime=@$ENV{SOURCE_DATE_EPOCH}" --clamp-mtime)
    endif()
    file(MAKE_DIRECTORY "${_work_dir}/member")
    execute_process(
      COMMAND "${_tar}" --create --file "${_work_dir}/member/control.tar" --format=gnu
              --owner=0 --group=0 --numeric-owner ${_mtime} --no-recursion
              -C "${_work_dir}/control" ${_control_members}
      RESULT_VARIABLE _result
      ERROR_VARIABLE _error
    )
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "Could not rebuild the control archive of ${_package}: ${_error}")
    endif()
    execute_process(
      COMMAND "${_gzip}" -9 --no-name "${_work_dir}/member/control.tar"
      RESULT_VARIABLE _result
      ERROR_VARIABLE _error
    )
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "Could not compress the control archive of ${_package}: ${_error}")
    endif()

    execute_process(
      COMMAND "${_ar}" t "${_package}"
      OUTPUT_VARIABLE _members
      RESULT_VARIABLE _result
      ERROR_VARIABLE _error
    )
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "Could not list the members of ${_package}: ${_error}")
    endif()
    if(NOT _members MATCHES "(^|\n)control\\.tar\\.gz\n")
      message(FATAL_ERROR "${_package} has no control.tar.gz member; its members are:\n${_members}")
    endif()

    ## 'ar r' replaces the member in place, keeping the order dpkg requires. ar writes the
    ## new archive to a temporary file of its own and renames it over the package, so an
    ## interrupted run leaves the package as it was. D: no time stamps or owners; S: no
    ## symbol table, which would become the first member.
    execute_process(
      COMMAND "${_ar}" rDS "${_package}" control.tar.gz
      WORKING_DIRECTORY "${_work_dir}/member"
      RESULT_VARIABLE _result
      ERROR_VARIABLE _error
    )
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "Could not replace the control archive of ${_package}: ${_error}")
    endif()

    ## Check the result as dpkg reads it.
    execute_process(
      COMMAND "${_ar}" t "${_package}"
      OUTPUT_VARIABLE _new_members
      RESULT_VARIABLE _result
    )
    execute_process(
      COMMAND "${_dpkg_deb}" --field "${_package}" Depends
      OUTPUT_VARIABLE _written_depends
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _field_result
      ERROR_VARIABLE _error
    )
    if(NOT _result EQUAL 0 OR NOT _new_members STREQUAL _members
       OR NOT _field_result EQUAL 0 OR NOT _written_depends STREQUAL _new_depends)
      message(FATAL_ERROR "The rewritten ${_package} does not read back as expected.\n"
                          "Members: ${_new_members}\nDepends: ${_written_depends}\n${_error}")
    endif()
    message(STATUS "cpack_deb_t64_dependencies.cmake: ${_package_name}: Depends: ${_new_depends}")
  endforeach()

  file(REMOVE_RECURSE "${_work_dir}")
endfunction()

_openms_deb_add_t64_alternatives("${CPACK_PACKAGE_DIRECTORY}/_openms_deb_t64" ${CPACK_PACKAGE_FILES})

cmake_policy(POP)

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------

## CPack runs this after writing each package and before copying it back to the
## build tree. dpkg-shlibdeps reports only the package names known to the build
## host, so replace Qt6 core/gui/widgets dependencies with alternatives for Ubuntu
## releases on either side of the t64 transition. Only the control archive is rebuilt;
## the compressed data archive is preserved byte-for-byte.
cmake_policy(VERSION 3.24)

if(NOT CPACK_PACKAGE_FILES)
  return()
endif()

find_program(_openms_dpkg_deb_EXECUTABLE dpkg-deb REQUIRED)
find_program(_openms_ar_EXECUTABLE ar REQUIRED)
find_program(_openms_tar_EXECUTABLE tar REQUIRED)

foreach(_openms_package IN LISTS CPACK_PACKAGE_FILES)
  if(NOT _openms_package MATCHES "\\.deb$")
    continue()
  endif()

  set(_openms_work_dir "${CPACK_PACKAGE_DIRECTORY}/_openms_deb_t64")
  file(REMOVE_RECURSE "${_openms_work_dir}")
  file(MAKE_DIRECTORY "${_openms_work_dir}/control")

  execute_process(
    COMMAND "${_openms_dpkg_deb_EXECUTABLE}" --ctrl-tarfile "${_openms_package}"
    OUTPUT_FILE "${_openms_work_dir}/control.tar"
    RESULT_VARIABLE _openms_control_result
    ERROR_VARIABLE _openms_control_error
  )
  if(NOT _openms_control_result EQUAL 0)
    message(FATAL_ERROR "Could not read control archive from ${_openms_package}: ${_openms_control_error}")
  endif()

  execute_process(
    COMMAND "${_openms_tar_EXECUTABLE}" -xf "${_openms_work_dir}/control.tar" -C "${_openms_work_dir}/control"
    RESULT_VARIABLE _openms_tar_result
    ERROR_VARIABLE _openms_tar_error
  )
  if(NOT _openms_tar_result EQUAL 0)
    message(FATAL_ERROR "Could not extract control archive from ${_openms_package}: ${_openms_tar_error}")
  endif()

  set(_openms_control_file "${_openms_work_dir}/control/control")
  file(READ "${_openms_control_file}" _openms_control)
  string(REGEX MATCH "(^|\n)Depends:[ \t]*([^\r\n]*)" _openms_depends_field "${_openms_control}")
  if(NOT _openms_depends_field)
    file(REMOVE_RECURSE "${_openms_work_dir}")
    continue()
  endif()

  set(_openms_depends_prefix "${CMAKE_MATCH_1}")
  set(_openms_depends "${CMAKE_MATCH_2}")
  string(REPLACE "," ";" _openms_depends_list "${_openms_depends}")
  set(_openms_rewritten_depends)
  foreach(_openms_dependency IN LISTS _openms_depends_list)
    string(STRIP "${_openms_dependency}" _openms_dependency)
    set(_openms_rewritten_dependency "${_openms_dependency}")
    foreach(_openms_base_package IN ITEMS libqt6core6 libqt6gui6 libqt6widgets6)
      if(_openms_dependency MATCHES "^${_openms_base_package}(t64)?([ \t]+\\([^)]*\\))?$")
        set(_openms_suffix "${CMAKE_MATCH_2}")
        if(CMAKE_MATCH_1)
          set(_openms_rewritten_dependency "${_openms_base_package}t64${_openms_suffix} | ${_openms_base_package}${_openms_suffix}")
        else()
          set(_openms_rewritten_dependency "${_openms_base_package}${_openms_suffix} | ${_openms_base_package}t64${_openms_suffix}")
        endif()
        break()
      endif()
    endforeach()
    list(APPEND _openms_rewritten_depends "${_openms_rewritten_dependency}")
  endforeach()
  list(JOIN _openms_rewritten_depends ", " _openms_rewritten_depends)

  if(NOT _openms_rewritten_depends STREQUAL _openms_depends)
    set(_openms_rewritten_field "${_openms_depends_prefix}Depends: ${_openms_rewritten_depends}")
    string(REPLACE "${_openms_depends_field}" "${_openms_rewritten_field}" _openms_control "${_openms_control}")
    file(WRITE "${_openms_control_file}" "${_openms_control}")

    execute_process(
      COMMAND "${_openms_tar_EXECUTABLE}" --owner=0 --group=0 --numeric-owner -cf "${_openms_work_dir}/control.tar" -C "${_openms_work_dir}/control" .
      RESULT_VARIABLE _openms_tar_result
      ERROR_VARIABLE _openms_tar_error
    )
    if(NOT _openms_tar_result EQUAL 0)
      message(FATAL_ERROR "Could not rebuild control archive for ${_openms_package}: ${_openms_tar_error}")
    endif()

    execute_process(
      COMMAND "${_openms_ar_EXECUTABLE}" t "${_openms_package}"
      OUTPUT_VARIABLE _openms_archive_members
      RESULT_VARIABLE _openms_ar_result
      ERROR_VARIABLE _openms_ar_error
    )
    if(NOT _openms_ar_result EQUAL 0)
      message(FATAL_ERROR "Could not inspect ${_openms_package}: ${_openms_ar_error}")
    endif()
    string(REPLACE "\r\n" "\n" _openms_archive_members "${_openms_archive_members}")
    string(REPLACE "\n" ";" _openms_archive_members "${_openms_archive_members}")
    set(_openms_data_member)
    foreach(_openms_member IN LISTS _openms_archive_members)
      if(_openms_member MATCHES "^data\\.tar(\\..+)?$")
        set(_openms_data_member "${_openms_member}")
        break()
      endif()
    endforeach()
    if(NOT _openms_data_member)
      message(FATAL_ERROR "Could not find data archive in ${_openms_package}")
    endif()

    execute_process(
      COMMAND "${_openms_ar_EXECUTABLE}" p "${_openms_package}" debian-binary
      OUTPUT_FILE "${_openms_work_dir}/debian-binary"
      RESULT_VARIABLE _openms_ar_result
      ERROR_VARIABLE _openms_ar_error
    )
    if(NOT _openms_ar_result EQUAL 0)
      message(FATAL_ERROR "Could not read package version from ${_openms_package}: ${_openms_ar_error}")
    endif()
    execute_process(
      COMMAND "${_openms_ar_EXECUTABLE}" p "${_openms_package}" "${_openms_data_member}"
      OUTPUT_FILE "${_openms_work_dir}/${_openms_data_member}"
      RESULT_VARIABLE _openms_ar_result
      ERROR_VARIABLE _openms_ar_error
    )
    if(NOT _openms_ar_result EQUAL 0)
      message(FATAL_ERROR "Could not read data archive from ${_openms_package}: ${_openms_ar_error}")
    endif()

    set(_openms_rebuilt_package "${_openms_package}.openms-t64")
    execute_process(
      COMMAND "${_openms_ar_EXECUTABLE}" qc "${_openms_rebuilt_package}"
        "${_openms_work_dir}/debian-binary"
        "${_openms_work_dir}/control.tar"
        "${_openms_work_dir}/${_openms_data_member}"
      RESULT_VARIABLE _openms_ar_result
      ERROR_VARIABLE _openms_ar_error
    )
    if(NOT _openms_ar_result EQUAL 0)
      message(FATAL_ERROR "Could not rebuild ${_openms_package}: ${_openms_ar_error}")
    endif()
    file(RENAME "${_openms_rebuilt_package}" "${_openms_package}" RESULT _openms_rename_result)
    if(NOT _openms_rename_result STREQUAL "0")
      message(FATAL_ERROR "Could not replace ${_openms_package}: ${_openms_rename_result}")
    endif()
  endif()

  file(REMOVE_RECURSE "${_openms_work_dir}")
endforeach()

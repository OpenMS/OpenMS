# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

## share/OpenMS/THIRD-PARTY-NOTICES.txt holds the license, copyright and notice texts of the
## third-party software an installation contains, in one file. It first lists the components,
## grouped by where they come from, each with the numbers of the texts that apply to it. The
## texts follow, each only once: components whose texts are identical (the same license and
## the same copyright holders) share a number. The files the texts come from stay where
## they are installed, in share/OpenMS/LICENSES and in the tool folders under
## share/OpenMS/THIRDPARTY.
##
## The code that installs these files records them with openms_add_third_party_notice():
## cmake/third_party_licenses.cmake, the Thermo bridge in cmake_findExternalLibs.cmake and
## the THIRDPARTY tools in package_general.cmake. openms_install_third_party_notices(), which
## the top-level CMakeLists.txt calls last, writes the list of recorded files when
## configuring, and the notices from that list when installing, in the component 'share'.
## The notices cannot be written when configuring, since some of the files (the Thermo
## bridge's) exist only after the build. They cannot be made from the installed tree either,
## because the package generators install each component into a folder of its own. The DEB
## also installs the notices, after the OpenMS license, as /usr/share/doc/openms/copyright,
## the file Debian's policy asks every package to have.
##
## The wheel builds add the licenses of the libraries that the wheel repair tools bundle
## after installing OpenMS (tools/ci/collect_wheel_licenses.py). They then run this file as a
## script, which writes the notices again, from the installed tree:
##   cmake -DSHARE_DIR=<prefix>/share/OpenMS -P cmake/third_party_notices.cmake

## The functions below keep these policies wherever they are called from: in script mode,
## and in the install script, which sets none. include() gives this file a policy scope of
## its own, so the setting does not leak into the including project.
cmake_policy(VERSION 3.24)

## openms_add_third_party_notice(<path> <file>)
##   Records that the installation contains <file> as share/OpenMS/<path>.
function(openms_add_third_party_notice path file)
  set_property(GLOBAL APPEND PROPERTY OPENMS_THIRD_PARTY_NOTICES "${path}|${file}")
endfunction()

## openms_add_third_party_notice_directory(<path> <directory>)
##   Records every file under <directory>, which is installed as share/OpenMS/<path>.
function(openms_add_third_party_notice_directory path directory)
  file(GLOB_RECURSE _files LIST_DIRECTORIES false RELATIVE "${directory}" "${directory}/*")
  foreach(_file IN LISTS _files)
    openms_add_third_party_notice("${path}/${_file}" "${directory}/${_file}")
  endforeach()
endfunction()

## Sets <var> to the license, copyright and notice files at the top of a tool folder, such as
## LICENSE, License.txt, NOTICE, COPYING or THIRD-PARTY-NOTICES.txt; READMEs are left out.
function(openms_third_party_notice_files_of_tool var directory)
  file(GLOB _files LIST_DIRECTORIES false RELATIVE "${directory}" "${directory}/*")
  set(_notices)
  foreach(_file IN LISTS _files)
    string(TOUPPER "${_file}" _name)
    if(_name MATCHES "(^|[-_])(LICEN[CS]E|COPYING|COPYRIGHT|NOTICES?|EULA)")
      list(APPEND _notices "${_file}")
    endif()
  endforeach()
  set(${var} "${_notices}" PARENT_SCOPE)
endfunction()

## openms_install_third_party_notices()
##   Installs share/OpenMS/THIRD-PARTY-NOTICES.txt, made from the files recorded so far.
function(openms_install_third_party_notices)
  set(_dir "${PROJECT_BINARY_DIR}/third_party_notices")
  get_property(_notices GLOBAL PROPERTY OPENMS_THIRD_PARTY_NOTICES)
  list(JOIN _notices "\n" _list)
  file(WRITE "${_dir}/files.txt" "${_list}\n")
  set(_debian_copyright "")
  ## PACKAGE_TYPE, not CPACK_GENERATOR: include(CPack) resets the latter to its defaults.
  if("${PACKAGE_TYPE}" STREQUAL "deb")
    set(_debian_copyright "LICENSE \"${OPENMS_HOST_DIRECTORY}/License.txt\" LICENSE_OUTPUT \"${_dir}/copyright\"")
  endif()
  install(CODE "
    include(\"${CMAKE_CURRENT_FUNCTION_LIST_FILE}\")
    openms_write_third_party_notices(\"${_dir}/THIRD-PARTY-NOTICES.txt\"
                                     FILE_LIST \"${_dir}/files.txt\" ${_debian_copyright})"
    COMPONENT share)
  install(FILES "${_dir}/THIRD-PARTY-NOTICES.txt"
          DESTINATION "${INSTALL_SHARE_DIR}"
          COMPONENT share)
  if(_debian_copyright)
    ## Debian policy 12.5: /usr/share/doc/<package>/copyright, neither compressed nor a
    ## symbolic link.
    string(TOLOWER "${CPACK_PACKAGE_NAME}" _package)
    install(FILES "${_dir}/copyright"
            DESTINATION "share/doc/${_package}"
            COMPONENT share)
  endif()
  list(LENGTH _notices _count)
  message(STATUS "THIRD-PARTY-NOTICES.txt: ${_count} license, copyright and notice files")
endfunction()

## Sets <text_var> to the content of <file> with Unix line ends and without trailing
## whitespace, and <binary_var> to TRUE if <file> is not a text file, i.e. has a NUL byte in
## its first 8 KiB (the Thermo license is a Word document). file(READ) drops carriage
## returns, so the length of what it reads says nothing.
function(_openms_third_party_notice_text text_var binary_var file)
  file(READ "${file}" _head LIMIT 8192 HEX)
  ## " <byte> <byte> ...", so that " 00" can only match a whole byte.
  string(REGEX REPLACE "(..)" " \\1" _head "${_head}")
  string(FIND "${_head}" " 00" _nul)
  if(NOT _nul EQUAL -1)
    set(${binary_var} TRUE PARENT_SCOPE)
    set(${text_var} "" PARENT_SCOPE)
    return()
  endif()
  file(READ "${file}" _text)
  string(REPLACE "\r" "" _text "${_text}")
  string(REGEX REPLACE "[ \t\n]+$" "" _text "${_text}")
  set(${binary_var} FALSE PARENT_SCOPE)
  set(${text_var} "${_text}" PARENT_SCOPE)
endfunction()

## Sets <group_var> to the number of the group <path> belongs to (the order of the groups in
## the notices) and <component_var> to the name of its component.
function(_openms_third_party_notice_component group_var component_var path)
  get_filename_component(_name "${path}" NAME)
  get_filename_component(_stem "${path}" NAME_WLE)
  if(path MATCHES "^LICENSES/Qt/")
    set(_group 1)
    set(_component "Qt")
  elseif(path MATCHES "^LICENSES/(vcpkg|homebrew|contrib|system|vendored)/(.+)$")
    set(_folder "${CMAKE_MATCH_1}")
    set(_rest "${CMAKE_MATCH_2}")
    set(_folders vcpkg homebrew contrib system vendored)
    list(FIND _folders "${_folder}" _index)
    math(EXPR _group "${_index} + 2")
    get_filename_component(_component "${_rest}" DIRECTORY)
    if(_name STREQUAL "INDEX.txt" AND NOT _component)
      set(_component "(the libraries each package provides)")
    elseif(NOT _component)
      ## vcpkg/<port>.txt
      set(_component "${_stem}")
    endif()
  elseif(path MATCHES "^LICENSES/ThermoRawFileReader|^openms_thermo_bridge/")
    set(_group 7)
    if(path MATCHES "^LICENSES/")
      set(_component "Thermo RawFileReader")
    else()
      set(_component "openms-thermo-bridge")
    endif()
  elseif(path MATCHES "^THIRDPARTY/([^/]+)/")
    set(_component "${CMAKE_MATCH_1}")
    set(_group 8)
  else()
    set(_group 9)
    set(_component "${path}")
  endif()
  set(${group_var} "${_group}" PARENT_SCOPE)
  set(${component_var} "${_component}" PARENT_SCOPE)
endfunction()

## openms_write_third_party_notices(<output> [FILE_LIST <list>] [SHARE_DIR <dir>]
##                                  [LICENSE <license> LICENSE_OUTPUT <output>])
##   Writes the notices for the files named in <list> (lines "<path>|<file>", as
##   openms_install_third_party_notices() writes them) and for those installed in the
##   share/OpenMS folder <dir>. With LICENSE, also writes LICENSE_OUTPUT: <license> followed
##   by the notices.
function(openms_write_third_party_notices output)
  cmake_parse_arguments(PARSE_ARGV 1 arg "" "FILE_LIST;SHARE_DIR;LICENSE;LICENSE_OUTPUT" "")
  set(_entries)
  if(arg_FILE_LIST)
    file(STRINGS "${arg_FILE_LIST}" _entries ENCODING UTF-8)
  endif()
  if(arg_SHARE_DIR)
    file(GLOB_RECURSE _files LIST_DIRECTORIES false RELATIVE "${arg_SHARE_DIR}"
         "${arg_SHARE_DIR}/LICENSES/*"
         "${arg_SHARE_DIR}/openms_thermo_bridge/managed/THIRD-PARTY-NOTICES.txt")
    ## OpenMS's own license is not third-party.
    list(REMOVE_ITEM _files "LICENSES/OpenMS-BSD-3-Clause.txt")
    file(GLOB _tools LIST_DIRECTORIES true RELATIVE "${arg_SHARE_DIR}/THIRDPARTY"
         "${arg_SHARE_DIR}/THIRDPARTY/*")
    foreach(_tool IN LISTS _tools)
      openms_third_party_notice_files_of_tool(_tool_files "${arg_SHARE_DIR}/THIRDPARTY/${_tool}")
      foreach(_file IN LISTS _tool_files)
        list(APPEND _files "THIRDPARTY/${_tool}/${_file}")
      endforeach()
    endforeach()
    foreach(_file IN LISTS _files)
      list(APPEND _entries "${_file}|${arg_SHARE_DIR}/${_file}")
    endforeach()
  endif()

  ## Order: by group, then by path; "<group>|<path>|<file>".
  set(_sorted)
  foreach(_entry IN LISTS _entries)
    string(REPLACE "|" ";" _entry "${_entry}")
    list(GET _entry 0 _path)
    list(GET _entry 1 _file)
    if(NOT EXISTS "${_file}")
      message(FATAL_ERROR "THIRD-PARTY-NOTICES.txt: ${_file} (installed as ${_path}) does not exist")
    endif()
    _openms_third_party_notice_component(_group _component "${_path}")
    list(APPEND _sorted "${_group}|${_path}|${_file}")
  endforeach()
  list(REMOVE_DUPLICATES _sorted)
  list(SORT _sorted)

  ## Numbers the distinct texts and collects, per component, the numbers of its texts.
  set(_hashes)        ## hash of each numbered text; its number is its index + 1
  set(_text_files)    ## the file each numbered text is read from
  set(_text_paths)    ## and where it is installed
  set(_components)    ## "<group>|<component>", in order
  foreach(_entry IN LISTS _sorted)
    string(REPLACE "|" ";" _entry "${_entry}")
    list(GET _entry 1 _path)
    list(GET _entry 2 _file)
    _openms_third_party_notice_component(_group _component "${_path}")
    set(_key "${_group}|${_component}")
    string(MD5 _refs "${_key}")
    set(_refs "_refs_${_refs}")
    if(NOT _key IN_LIST _components)
      list(APPEND _components "${_key}")
      set(${_refs} "")
    endif()
    _openms_third_party_notice_text(_text _binary "${_file}")
    if(_binary)
      string(APPEND ${_refs} " ${_path} (not a text file)")
      continue()
    endif()
    ## Texts that differ only in white space (line breaks, indentation) are the same text.
    string(REGEX REPLACE "[ \t\n]+" " " _hash "${_text}")
    string(SHA256 _hash "${_hash}")
    list(FIND _hashes "${_hash}" _index)
    if(_index EQUAL -1)
      list(APPEND _hashes "${_hash}")
      list(APPEND _text_files "${_file}")
      list(APPEND _text_paths "${_path}")
      list(LENGTH _hashes _number)
    else()
      math(EXPR _number "${_index} + 1")
    endif()
    if(NOT "${${_refs}}" MATCHES " \\[${_number}\\]( |$)")
      string(APPEND ${_refs} " [${_number}]")
    endif()
  endforeach()

  set(_titles
    "Qt"
    "Libraries from vcpkg"
    "Libraries from Homebrew"
    "Libraries from the OpenMS contrib"
    "Libraries of the build system"
    "Third-party code compiled into OpenMS"
    "Thermo RawFileReader"
    "Tools in share/OpenMS/THIRDPARTY"
    "Other")
  set(_rule "==========================================================================================")
  string(CONCAT _header
    "THIRD-PARTY SOFTWARE NOTICES\n"
    "\n"
    "This installation of OpenMS contains the third-party software listed below. Each component\n"
    "is followed by the numbers of the license, copyright and notice texts that apply to it. The\n"
    "texts come after the list, each only once: components whose texts are identical share a\n"
    "number. Paths are relative to share/OpenMS, the folder in which the installation keeps the\n"
    "same files as their projects distribute them: LICENSES, and the folders of the tools under\n"
    "THIRDPARTY. In a pyOpenMS wheel, that folder is pyopenms/share/OpenMS. The license of\n"
    "OpenMS itself is LICENSES/OpenMS-BSD-3-Clause.txt.\n")
  set(_notices "${output}.tmp")
  file(WRITE "${_notices}" "${_header}")
  set(_current_group "")
  foreach(_key IN LISTS _components)
    string(REPLACE "|" ";" _parts "${_key}")
    list(GET _parts 0 _group)
    list(GET _parts 1 _component)
    if(NOT _group STREQUAL _current_group)
      math(EXPR _title_index "${_group} - 1")
      list(GET _titles ${_title_index} _title)
      file(APPEND "${_notices}" "\n${_title}\n")
      set(_current_group "${_group}")
    endif()
    string(MD5 _refs "${_key}")
    set(_refs "_refs_${_refs}")
    string(LENGTH "${_component}" _length)
    set(_padding " ")
    if(_length LESS 44)
      math(EXPR _pad "44 - ${_length}")
      string(REPEAT " " ${_pad} _padding)
    endif()
    file(APPEND "${_notices}" "  ${_component}${_padding}${${_refs}}\n")
  endforeach()

  set(_number 0)
  foreach(_file IN LISTS _text_files)
    list(GET _text_paths ${_number} _path)
    math(EXPR _number "${_number} + 1")
    _openms_third_party_notice_text(_text _binary "${_file}")
    file(APPEND "${_notices}" "\n\n${_rule}\n[${_number}] ${_path}\n${_rule}\n\n${_text}\n")
  endforeach()
  file(RENAME "${_notices}" "${output}")

  if(arg_LICENSE)
    file(READ "${arg_LICENSE}" _license)
    file(READ "${output}" _notices_text)
    file(WRITE "${arg_LICENSE_OUTPUT}" "${_license}\n\n${_rule}\n\n${_notices_text}")
  endif()

  list(LENGTH _components _component_count)
  list(LENGTH _text_files _text_count)
  message(STATUS "Wrote ${output}: ${_component_count} components, ${_text_count} texts")
endfunction()

## Script mode, for the wheels: cmake -DSHARE_DIR=<prefix>/share/OpenMS -P <this file>
if(CMAKE_SCRIPT_MODE_FILE STREQUAL CMAKE_CURRENT_LIST_FILE)
  if(NOT SHARE_DIR OR NOT IS_DIRECTORY "${SHARE_DIR}")
    message(FATAL_ERROR "Pass the installed share/OpenMS folder as -DSHARE_DIR=<folder>.")
  endif()
  openms_write_third_party_notices("${SHARE_DIR}/THIRD-PARTY-NOTICES.txt" SHARE_DIR "${SHARE_DIR}")
endif()

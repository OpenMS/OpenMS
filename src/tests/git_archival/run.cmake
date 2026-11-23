# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.24)

# cmake -DBINARY_DIR=<scratch> -P run.cmake
#
# Tests git_archival_info() (cmake/Modules/GetGitRevisionDescription.cmake), from which a source
# tree made by git archive, such as the release tarball, takes its version: on .git_archival.txt
# as a checkout has it and as git archive fills it in. The cases do not read the real file, which
# is filled in when this runs from a source tarball.
if(NOT BINARY_DIR)
  message(FATAL_ERROR "Set BINARY_DIR to a scratch directory")
endif()
file(MAKE_DIRECTORY "${BINARY_DIR}")
include("${CMAKE_CURRENT_LIST_DIR}/../../../cmake/Modules/GetGitRevisionDescription.cmake")

# _check(<case> <file content> <found> <describe> <hash> <date>)
function(_check case content exp_found exp_describe exp_hash exp_date)
  set(file "${BINARY_DIR}/${case}.txt")
  file(WRITE "${file}" "${content}")
  set(describe "<unset>")
  set(hash "<unset>")
  set(date "<unset>")
  git_archival_info("${file}" describe hash date found)
  set(got "${found}|${describe}|${hash}|${date}")
  set(want "${exp_found}|${exp_describe}|${exp_hash}|${exp_date}")
  if(got STREQUAL want)
    message(STATUS "${case}: ok")
  else()
    message(SEND_ERROR "${case}: got '${got}', expected '${want}'")
  endif()
endfunction()

set(_node "5d5cbff4053b281763a1a79bf69e81c27967cfdf")
set(_placeholder "%(describe:tags=true,abbrev=40,match=v[0-9]*,exclude=v*[!0-9.]*)")

# a checkout: the placeholders are not filled in
_check(checkout "node: $Format:%H$\nnode-date: $Format:%cI$\ndescribe-name: $Format:${_placeholder}$\n"
  FALSE "<unset>" "<unset>" "<unset>")
# the archive of a release tag
_check(tag "node: ${_node}\nnode-date: 2026-09-30T08:25:00+02:00\ndescribe-name: v3.6.0\n"
  TRUE "v3.6.0" "5d5cbff" "2026-09-30 08:25:00 +0200")
# a UTC commit date: git 2.45 and later write Z, older git +00:00
_check(utc_z "node: ${_node}\nnode-date: 2026-09-30T06:25:00Z\ndescribe-name: v3.6.0\n"
  TRUE "v3.6.0" "5d5cbff" "2026-09-30 06:25:00 +0000")
_check(utc_offset "node: ${_node}\nnode-date: 2026-09-30T06:25:00+00:00\ndescribe-name: v3.6.0\n"
  TRUE "v3.6.0" "5d5cbff" "2026-09-30 06:25:00 +0000")
# an untagged commit, and one without a release tag in its history (e.g. from a depth-1 clone)
_check(untagged "node: ${_node}\nnode-date: 2026-09-30T06:25:00Z\ndescribe-name: v3.6.0-6-g${_node}\n"
  TRUE "v3.6.0-6-g${_node}" "5d5cbff" "2026-09-30 06:25:00 +0000")
_check(no_tag "node: ${_node}\nnode-date: 2026-09-30T06:25:00Z\ndescribe-name: \n"
  TRUE "" "5d5cbff" "2026-09-30 06:25:00 +0000")
# git before 2.35 fills in node, but leaves the describe placeholder as it is (with a warning)
_check(old_git "node: ${_node}\nnode-date: 2026-09-30T06:25:00+00:00\ndescribe-name: ${_placeholder}\n"
  TRUE "" "5d5cbff" "2026-09-30 06:25:00 +0000")

# no .git_archival.txt at all
set(found "<unset>")
git_archival_info("${BINARY_DIR}/does-not-exist.txt" describe hash date found)
if(found STREQUAL "FALSE")
  message(STATUS "missing: ok")
else()
  message(SEND_ERROR "missing: got '${found}', expected 'FALSE'")
endif()

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# $Maintainer: Timo Sachsenberg $

foreach(_argument IN ITEMS EXECUTE_PIPELINE WORKFLOW RESOURCE_FILE TEST_DIRECTORY)
  if(NOT DEFINED ${_argument} OR "${${_argument}}" STREQUAL "")
    message(FATAL_ERROR "Missing required test argument ${_argument}")
  endif()
endforeach()

# A fresh child directory prevents stale results from satisfying assertions on
# subsequent CTest runs, and preserves failed-run artifacts for investigation.
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_directory "${TEST_DIRECTORY}/run_${_run_id}")
set(_working_directory "${_run_directory}/work")
set(_home_directory "${_run_directory}/home")
file(MAKE_DIRECTORY "${_working_directory}" "${_home_directory}")
file(COPY_FILE "${WORKFLOW}" "${_working_directory}/analysis.v2.toppas")
file(COPY_FILE "${RESOURCE_FILE}" "${_working_directory}/inputs.trf")

# The workflow's embedded input lists are empty; a successful run therefore
# verifies that the configured resource file supplies all input branches.
# Deliberately omit -out_dir and use a filename containing multiple dots.
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env "OPENMS_HOME_PATH=${_home_directory}"
    "${EXECUTE_PIPELINE}" -test -in analysis.v2.toppas -resource_file inputs.trf -num_jobs 2
  WORKING_DIRECTORY "${_working_directory}"
  TIMEOUT 150
  RESULT_VARIABLE _exit_code
  OUTPUT_VARIABLE _stdout
  ERROR_VARIABLE _stderr)
if(NOT "${_exit_code}" STREQUAL "0")
  message(FATAL_ERROR "ExecutePipeline failed (${_exit_code}).\n${_stdout}\n${_stderr}")
endif()

# Historical QFileInfo::baseName() chooses the portion before the first dot.
set(_output_directory "${_home_directory}/analysis")
if(NOT EXISTS "${_output_directory}/TOPPAS.log")
  message(FATAL_ERROR "Missing default-output log: ${_output_directory}/TOPPAS.log\n${_stdout}\n${_stderr}")
endif()
if(EXISTS "${_home_directory}/analysis.v2")
  message(FATAL_ERROR "Default output directory used the last-dot stem instead of the historical first-dot basename.")
endif()

# Two FileInfo branches process two files and one file respectively; the
# collector feeds both files to a single FileMerger invocation.
set(_output_subdirectories "007-FileInfo-out_tsv" "008-FileInfo-out_tsv" "010-FileMerger-out")
set(_expected_counts 2 1 1)
foreach(_index RANGE 0 2)
  list(GET _output_subdirectories ${_index} _subdirectory)
  list(GET _expected_counts ${_index} _expected_count)
  set(_directory "${_output_directory}/TOPPAS_out/${_subdirectory}")
  file(GLOB _files LIST_DIRECTORIES FALSE "${_directory}/*")
  list(LENGTH _files _count)
  if(NOT _count EQUAL _expected_count)
    message(FATAL_ERROR "Expected ${_expected_count} published files in ${_directory}, found ${_count}.\n${_stdout}\n${_stderr}")
  endif()
  foreach(_file IN LISTS _files)
    file(SIZE "${_file}" _size)
    if(_size EQUAL 0)
      message(FATAL_ERROR "Published result is empty: ${_file}")
    endif()
  endforeach()
endforeach()

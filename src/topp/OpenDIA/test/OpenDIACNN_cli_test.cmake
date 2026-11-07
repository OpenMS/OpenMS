# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

file(MAKE_DIRECTORY "${WORK}/input.d")
file(COPY_FILE "${LIBRARY}" "${WORK}/library.pqp")
file(WRITE "${WORK}/input.d/analysis.tdf" "raw input sentinel\n")
file(WRITE "${WORK}/windows.tsv" "sample_id\tprecursor_id\tgroup_id\tlabel\trt_apex\n")
foreach(model IN ITEMS "${WORK}/library.pqp" "${WORK}/input.d/analysis.tdf")
  file(SHA256 "${model}" before)
  execute_process(COMMAND "${BINARY}" -mode train -in "${WORK}/input.d"
    -lib "${WORK}/library.pqp" -windows "${WORK}/windows.tsv"
    -model "${model}" -out "${WORK}/scores.tsv"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "must not overwrite")
    message(FATAL_ERROR "Native CLI did not reject input overwrite: ${status}\n${output}\n${error}")
  endif()
  file(SHA256 "${model}" after)
  if(NOT before STREQUAL after OR EXISTS "${WORK}/scores.tsv")
    message(FATAL_ERROR "Native CLI changed input or opened score output before validation")
  endif()
endforeach()

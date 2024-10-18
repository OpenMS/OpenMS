# Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# include helper functions 
include ( ${SCRIPT_DIR}common.cmake )

set(required_variables "TOOLNAME;PARAM;CTD_PATH")
check_variables(required_variables)

set(filename "${CTD_PATH}/${TOOLNAME}.ctd")

file(READ "${filename}" contents)
string(REPLACE "[" "_OPENBRACKET_" contents "${contents}")
string(REPLACE "]" "_CLOSEBRACKET_" contents "${contents}")
string(REPLACE ";" "\\;" contents "${contents}")
string(REGEX REPLACE "\n" ";" contents "${contents}")

set(APPEND_TO_FILE FALSE)
foreach(line ${contents})
  set(pos -1) 
  if(line MATCHES ".*name=\"${PARAM}\".*") 
    string(LENGTH "${CMAKE_MATCH_1}" pos) 
  endif()
  
  string(REPLACE "_OPENBRACKET_" "[" line "${line}")
  string(REPLACE "_CLOSEBRACKET_" "]" line "${line}")
  
  # we only write out line that do not contain our parameter
  if (pos EQUAL -1)
    if(APPEND_TO_FILE)
      file(APPEND "${filename}" "${line}\n")
    else() 
      # write the first line
      file(WRITE "${filename}" "${line}\n")
      # append the rest
      set(APPEND_TO_FILE TRUE)
    endif()
  endif()
endforeach()

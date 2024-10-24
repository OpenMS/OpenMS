# Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# include helper functions 
include ( ${SCRIPT_DIR}common.cmake )

set(required_variables "TOOLNAME;CTD_FILE")
check_variables(required_variables)

set(filename ${CTD_FILE})

file(READ "${filename}" contents)
string(REPLACE "[" "_OPENBRACKET_" contents "${contents}")
string(REPLACE "]" "_CLOSEBRACKET_" contents "${contents}")
string(REPLACE ";" "\\\;" contents "${contents}")
string(REGEX REPLACE "\n" ";" contents "${contents}")

set(APPEND_TO_FILE FALSE)
foreach(line ${contents})
  string(REGEX REPLACE "(.*tool.*)name=\"([^ ]*)\"(.*)" "\\1executableName=\"\\2\" name=\"${TOOLNAME}\"\\3" line ${line})
  
  string(REPLACE "_OPENBRACKET_" "[" line "${line}")
  string(REPLACE "_CLOSEBRACKET_" "]" line "${line}")

  if(APPEND_TO_FILE)
    file(APPEND "${filename}" "${line}\n")
  else() 
    # write the first line
    file(WRITE "${filename}" "${line}\n")
    # append the rest
    set(APPEND_TO_FILE TRUE)
  endif()
endforeach()

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# include helper functions 
include ( ${SCRIPT_DIR}common.cmake )

set(required_variables "SOURCE_PATH;TARGET_DIRECTORY")
check_variables(required_variables)

file(GLOB_RECURSE share_files "${SOURCE_PATH}/share/*")

foreach(file ${share_files})  
  # remove the prefix 
  string(REPLACE "${SOURCE_PATH}/share/" "" trimmed_file ${file})
  
  set(pos -1) 
  if(trimmed_file MATCHES ".git/.*") 
    string(LENGTH "${CMAKE_MATCH_1}" pos) 
  endif() 
  
  # we only copy files that do not match
  if (pos EQUAL -1)  
    # copy
    configure_file(${file} ${TARGET_DIRECTORY}/${trimmed_file} COPYONLY)
  endif()
endforeach() 

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

set(required_variables "SE_PATH;TARGET_DIRECTORY")
check_variables(required_variables)

file(GLOB_RECURSE se_files "${SE_PATH}/*")
foreach(file ${se_files})  
  # remove the prefix 
  string(REPLACE "${SE_PATH}/" "" trimmed_file ${file})
  
  set(pos -1) 
  ## TODO currently we omit the vcredists and request the user to install
  ## them via the OpenMS dependency installer during plugin install in KNIME
  if(trimmed_file MATCHES ".git/.*" OR trimmed_file MATCHES "vcredist.*")
    string(LENGTH "${CMAKE_MATCH_1}" pos) 
  endif() 
  
  # we only write out line that do not contain our parameter
  if (pos EQUAL -1)  
    # copy
    configure_file(${file} ${TARGET_DIRECTORY}/${trimmed_file} COPYONLY)
  endif()
endforeach() 

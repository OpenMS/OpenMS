# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Timo Sachsenberg $

cmake_minimum_required(VERSION 3.24)

# Invoked by CXX_COMPILER_LAUNCHER. Bypass ccache only for /Yc; ccache can cache
# /Yu compilations and hashes the PCH input. This keeps machine-specific MSVC
# PCH artifacts out of shared compiler caches without disabling object caching.
set(_after_separator OFF)
set(_create_pch OFF)
math(EXPR _last "${CMAKE_ARGC} - 1")
foreach(_index RANGE 0 ${_last})
  set(_arg "${CMAKE_ARGV${_index}}")
  if(NOT _after_separator)
    if(_arg STREQUAL "--")
      set(_after_separator ON)
    endif()
    continue()
  endif()
  string(REPLACE ";" "\\;" _escaped "${_arg}")
  list(APPEND _command "${_escaped}")
  set(_options "${_arg}")
  if(_arg MATCHES "^@(.+)$")
    file(READ "${CMAKE_MATCH_1}" _response)
    separate_arguments(_options WINDOWS_COMMAND "${_response}")
  endif()
  foreach(_option IN LISTS _options)
    if(_option MATCHES "^[-/]Yc")
      set(_create_pch ON)
    endif()
  endforeach()
endforeach()

if(_create_pch)
  list(SUBLIST _command ${OPENMS_LAUNCHER_COUNT} -1 _command)
endif()
execute_process(COMMAND ${_command} RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "Compilation failed (${_result})")
endif()

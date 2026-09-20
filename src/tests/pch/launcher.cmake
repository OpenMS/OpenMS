# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Timo Sachsenberg $

# Mock compiler/cache used to test the MSVC launcher on every host. Record the
# command exactly, including spaces, semicolons and response-file arguments.
file(APPEND "${TRACE}" "${ROLE}\n")
set(_after_separator OFF)
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
  if(ROLE STREQUAL "compiler")
    file(APPEND "${TRACE}" "[${_arg}]\n")
  endif()
endforeach()
if(ROLE STREQUAL "cache")
  execute_process(COMMAND ${_command} RESULT_VARIABLE _result)
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "Mock compiler failed")
  endif()
endif()

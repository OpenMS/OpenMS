# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.15)

if(NOT DEFINED OUTPUT_DIR)
  message(FATAL_ERROR "OUTPUT_DIR must be defined")
endif()

if(NOT DEFINED TOPP_EXECUTABLES_CMAKE)
  message(FATAL_ERROR "TOPP_EXECUTABLES_CMAKE must be defined")
endif()

# Load the existing list of TOPP tools from executables.cmake
include("${TOPP_EXECUTABLES_CMAKE}")

# Combine all TOPP executables (regular and GUI-based)
set(_tool_names ${TOPP_executables} ${TOPP_executables_with_GUIlib})

if(_tool_names)
  list(REMOVE_DUPLICATES _tool_names)
  list(SORT _tool_names)
endif()

file(MAKE_DIRECTORY "${OUTPUT_DIR}")

foreach(_tool ${_tool_names})
  set(_cli_file "${OUTPUT_DIR}/TOPP_${_tool}.cli")
  if(NOT EXISTS "${_cli_file}")
    file(WRITE "${_cli_file}" "The TOPP tool '${_tool}' is not available in this build.\n")
  endif()

  set(_html_file "${OUTPUT_DIR}/TOPP_${_tool}.html")
  if(NOT EXISTS "${_html_file}")
    file(WRITE "${_html_file}"
         "<html><body><h2>Tool not available</h2><p>The TOPP tool '${_tool}' is not part of this build configuration.</p></body></html>\n")
  endif()
endforeach()

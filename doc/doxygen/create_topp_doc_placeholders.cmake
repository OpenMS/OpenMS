# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.15)

if(NOT DEFINED OUTPUT_DIR)
  message(FATAL_ERROR "OUTPUT_DIR must be defined")
endif()

if(NOT DEFINED SEARCH_ROOTS)
  message(FATAL_ERROR "SEARCH_ROOTS must be defined")
endif()

set(_tool_names)

foreach(_root ${SEARCH_ROOTS})
  if(NOT EXISTS "${_root}")
    message(WARNING "TOPP documentation placeholder search path '${_root}' does not exist.")
    continue()
  endif()

  file(GLOB_RECURSE _candidates LIST_DIRECTORIES FALSE
    "${_root}/*.cpp"
    "${_root}/*.h"
    "${_root}/*.hpp"
  )

  foreach(_candidate ${_candidates})
    file(STRINGS "${_candidate}" _doc_lines REGEX "@page[ \t]+TOPP_[A-Za-z0-9_]+")
    foreach(_line ${_doc_lines})
      string(REGEX MATCH "@page[ \t]+TOPP_[A-Za-z0-9_]+" _match "${_line}")
      if(_match)
        string(REGEX REPLACE "^@page[ \t]+TOPP_" "" _tool "${_match}")
        list(APPEND _tool_names "${_tool}")
      endif()
    endforeach()
  endforeach()
endforeach()

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

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Chris Bielow $
# $Authors: Chris Bielow $
# --------------------------------------------------------------------------

## Run with 'cmake -P', which starts with no policy settings of its own; IN_LIST below needs
## CMP0057. Matches the minimum the project requires.
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

### Checks the generated TOPP tool registry.
###
### The registry is written by openms_write_tool_registry() (cmake/topp_tool_macros.cmake) from
### the same openms_topp_tool() declarations that produce the executables, so the registry and
### the build cannot describe different sets of tools; that is the point of generating it and
### there is nothing to cross-check. What generating it does not rule out is the generator
### writing a file ToolHandler cannot read back, or losing entries on the way: either costs the
### tools in it their registration, and the only symptom is tools refusing to start.
###
### Driven by src/tests/topp/CMakeLists.txt. Arguments:
###   REGISTRY_FILE   the generated registry
###   TOOLS           the tools this build produces (${TOPP_TOOLS}), ';'-separated

foreach(_required REGISTRY_FILE TOOLS)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "check_tool_registry.cmake: -D${_required}= is required")
  endif()
endforeach()

if(NOT EXISTS "${REGISTRY_FILE}")
  message(FATAL_ERROR "The tool registry '${REGISTRY_FILE}' was not generated.")
endif()

file(STRINGS "${REGISTRY_FILE}" _lines)

## Each entry is '<name><TAB><category>', optionally followed by a TAB and the tool's '-type'
## sub-modes. ToolHandler reads it with CsvFile, which drops '#' comment lines; a line that is
## neither a comment nor a well-formed entry is one ToolHandler would log and skip, so the
## tool it was meant to register would silently not exist.
set(registered)
set(_line_number 0)
foreach(_line IN LISTS _lines)
  math(EXPR _line_number "${_line_number} + 1")
  if(_line MATCHES "^#" OR _line STREQUAL "")
    continue()
  endif()
  ## CMake cannot match a literal tab portably inside a bracket expression, so split on it.
  string(REPLACE "\t" ";" _fields "${_line}")
  list(LENGTH _fields _field_count)
  if(_field_count LESS 2)
    message(FATAL_ERROR "'${REGISTRY_FILE}' line ${_line_number} is not a registry entry: '${_line}'. "
                        "Expected '<tool name><TAB><category>'; ToolHandler skips such a line, so the tool "
                        "would not be registered at all.")
  endif()
  list(GET _fields 0 _name)
  list(GET _fields 1 _category)
  if(_name STREQUAL "" OR _category STREQUAL "")
    message(FATAL_ERROR "'${REGISTRY_FILE}' line ${_line_number} has an empty name or category: '${_line}'.")
  endif()
  if(_name IN_LIST registered)
    message(FATAL_ERROR "'${_name}' is registered twice in '${REGISTRY_FILE}'. ToolHandler throws on a duplicate "
                        "name, which stops every TOPP tool from starting.")
  endif()
  list(APPEND registered "${_name}")
endforeach()

## The generated file has to carry the declarations across intact. This cannot fail through a
## forgotten declaration (the same list builds the executables), but it does catch the
## generator dropping or truncating entries.
set(built ${TOOLS})
list(REMOVE_DUPLICATES built)

set(missing ${built})
list(REMOVE_ITEM missing ${registered})
if(missing)
  string(REPLACE ";" ", " _list "${missing}")
  message(FATAL_ERROR "The generated registry '${REGISTRY_FILE}' is missing tools this build declares: ${_list}. "
                      "openms_write_tool_registry() lost them.")
endif()

set(unexpected ${registered})
list(REMOVE_ITEM unexpected ${built})
if(unexpected)
  string(REPLACE ";" ", " _list "${unexpected}")
  message(FATAL_ERROR "The generated registry '${REGISTRY_FILE}' lists tools this build does not declare: ${_list}. "
                      "A registered tool without a binary is offered by TOPPAS and cannot run.")
endif()

list(LENGTH registered _count)
message(STATUS "Generated TOPP tool registry is readable and carries all ${_count} declared tools.")

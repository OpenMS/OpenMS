# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# --------------------------------------------------------------------------

include_guard(GLOBAL)

option(OPENMS_USE_PCH "Use private precompiled standard-library headers for OpenMS libraries" OFF)

function(openms_add_precompiled_headers target_name)
  if(NOT OPENMS_USE_PCH)
    return()
  endif()
  get_target_property(_pch_disabled ${target_name} DISABLE_PRECOMPILE_HEADERS)
  if(_pch_disabled)
    return()
  endif()

  # Build a separate PCH for each library and configuration: DLL export defines,
  # runtime libraries and compiler flags can differ. Do not propagate PCHs to
  # consumers or create a PCH for each of the many single-source tools/tests.
  # Keep this list independent of OpenMS and third-party configuration macros.
  foreach(_header IN ITEMS algorithm array cstddef cstdint functional iterator
      limits map memory numeric set stdexcept string string_view tuple
      type_traits unordered_map unordered_set utility vector)
    target_precompile_headers(${target_name} PRIVATE
      "$<$<COMPILE_LANGUAGE:CXX>:<${_header}$<ANGLE-R>>")
  endforeach()

  # Build MSVC PCHs locally (Microsoft does not guarantee cross-machine reuse),
  # while retaining ccache for the object files that consume them. Other
  # launchers are left alone, and no global cache configuration is changed.
  get_target_property(_launcher ${target_name} CXX_COMPILER_LAUNCHER)
  if(_launcher)
    list(GET _launcher 0 _launcher_command)
    get_filename_component(_launcher_name "${_launcher_command}" NAME_WE)
    if(_launcher_name STREQUAL "ccache")
      if(MSVC)
        execute_process(COMMAND ${_launcher} --version
          OUTPUT_VARIABLE _ccache_version RESULT_VARIABLE _ccache_result)
        if(NOT _ccache_result EQUAL 0 OR NOT _ccache_version MATCHES "ccache version ([0-9]+\\.[0-9]+(\\.[0-9]+)?)")
          message(FATAL_ERROR "Cannot determine the ccache version for MSVC PCH support")
        endif()
        if(CMAKE_MATCH_1 VERSION_LESS "4.14")
          message(FATAL_ERROR "MSVC PCH with ccache requires ccache 4.14 or newer. Upgrade ccache or configure without the C++ compiler launcher.")
        endif()
        list(LENGTH _launcher _launcher_count)
        set_property(TARGET ${target_name} PROPERTY CXX_COMPILER_LAUNCHER
          "${CMAKE_COMMAND};-DOPENMS_LAUNCHER_COUNT=${_launcher_count};-P;${CMAKE_CURRENT_FUNCTION_LIST_DIR}/msvc_pch_launcher.cmake;--;${_launcher}")
      elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${target_name} PRIVATE
          "$<$<COMPILE_LANGUAGE:CXX>:-Xclang;-fno-pch-timestamp>")
      endif()
    endif()
  endif()
endfunction()

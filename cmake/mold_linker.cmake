# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: OpenMS Development Team $
# $Authors: OpenMS Development Team $
# --------------------------------------------------------------------------

# This module detects and enables the mold linker if available and requested.
# mold is a modern, high-performance linker that can significantly speed up
# link times compared to traditional linkers like ld or gold.
#
# Usage:
#   Set USE_MOLD_LINKER to ON to enable mold linker detection and usage.
#
# Variables set by this module:
#   MOLD_FOUND - TRUE if mold linker was found
#   MOLD_LINKER_FLAGS - Linker flags to use mold

option(USE_MOLD_LINKER "Use mold linker if available (Linux only)" OFF)

if(USE_MOLD_LINKER)
  # mold is primarily a Linux linker
  if(UNIX AND NOT APPLE)
    # Define the linker flag as a constant
    set(_MOLD_LINKER_FLAG "-fuse-ld=mold")
    
    # Check if mold is available
    find_program(MOLD_EXECUTABLE NAMES mold)
    
    if(MOLD_EXECUTABLE)
      message(STATUS "Found mold linker: ${MOLD_EXECUTABLE}")
      
      # Check if the compiler supports using mold
      # -fuse-ld=mold is supported by GCC 9+ and Clang 9+
      include(CheckLinkerFlag)
      
      # Use CMAKE_REQUIRED_LINK_OPTIONS for linker flags (CMake 3.18+)
      # Fall back to check_cxx_compiler_flag for older CMake versions
      if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.18")
        check_linker_flag(CXX "${_MOLD_LINKER_FLAG}" COMPILER_SUPPORTS_MOLD)
      else()
        include(CheckCXXCompilerFlag)
        set(_ORIGINAL_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
        set(CMAKE_REQUIRED_FLAGS "${_MOLD_LINKER_FLAG}")
        check_cxx_compiler_flag("${_MOLD_LINKER_FLAG}" COMPILER_SUPPORTS_MOLD)
        set(CMAKE_REQUIRED_FLAGS "${_ORIGINAL_CMAKE_REQUIRED_FLAGS}")
        unset(_ORIGINAL_CMAKE_REQUIRED_FLAGS)
      endif()
      
      if(COMPILER_SUPPORTS_MOLD)
        message(STATUS "Enabling mold linker")
        set(MOLD_FOUND TRUE)
        set(MOLD_LINKER_FLAGS "${_MOLD_LINKER_FLAG}")
        
        # Apply to all linker types using add_link_options()
        # This is the modern CMake way (CMake 3.13+) and avoids duplication
        add_link_options("${MOLD_LINKER_FLAGS}")
      else()
        message(WARNING "Compiler does not support ${_MOLD_LINKER_FLAG} flag")
        set(MOLD_FOUND FALSE)
      endif()
    else()
      message(WARNING "mold linker requested but not found in PATH")
      set(MOLD_FOUND FALSE)
    endif()
    
    unset(_MOLD_LINKER_FLAG)
  else()
    message(STATUS "mold linker is only supported on Linux, skipping")
    set(MOLD_FOUND FALSE)
  endif()
endif()

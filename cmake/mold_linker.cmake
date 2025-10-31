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
    # Check if mold is available
    find_program(MOLD_EXECUTABLE NAMES mold)
    
    if(MOLD_EXECUTABLE)
      message(STATUS "Found mold linker: ${MOLD_EXECUTABLE}")
      
      # Check if the compiler supports using mold
      # We use -fuse-ld=mold which is supported by GCC 12.1+ and Clang 12+
      include(CheckCXXCompilerFlag)
      
      # Save original CMAKE_REQUIRED_FLAGS and restore it after the check
      set(_ORIGINAL_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
      set(CMAKE_REQUIRED_FLAGS "-fuse-ld=mold")
      check_cxx_compiler_flag("-fuse-ld=mold" COMPILER_SUPPORTS_MOLD)
      set(CMAKE_REQUIRED_FLAGS "${_ORIGINAL_CMAKE_REQUIRED_FLAGS}")
      unset(_ORIGINAL_CMAKE_REQUIRED_FLAGS)
      
      if(COMPILER_SUPPORTS_MOLD)
        message(STATUS "Enabling mold linker")
        set(MOLD_FOUND TRUE)
        set(MOLD_LINKER_FLAGS "-fuse-ld=mold")
        
        # Apply to all linker types
        add_link_options("${MOLD_LINKER_FLAGS}")
        
        # Also set the linker flags variables for older CMake compatibility
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MOLD_LINKER_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${MOLD_LINKER_FLAGS}")
        set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${MOLD_LINKER_FLAGS}")
      else()
        message(WARNING "Compiler does not support -fuse-ld=mold flag")
        set(MOLD_FOUND FALSE)
      endif()
    else()
      message(WARNING "mold linker requested but not found in PATH")
      set(MOLD_FOUND FALSE)
    endif()
  else()
    message(STATUS "mold linker is only supported on Linux, skipping")
    set(MOLD_FOUND FALSE)
  endif()
endif()

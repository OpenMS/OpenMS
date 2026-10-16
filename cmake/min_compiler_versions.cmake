# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Chris Bielow $
# $Authors: Chris Bielow $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Single source of truth for the minimum supported compiler versions.
#
# The values below are consumed in two places, so they only ever need to be
# edited here:
#   1. the checks at the bottom of this file, which abort configuration with a
#      readable message on an unsupported compiler
#   2. doc/doxygen/Doxyfile.in, which turns them into doxygen ALIASES
#      (@minMSVC, @minGCC, ...) so that the install documentation quotes these
#      exact numbers instead of hardcoding its own copy
#
# Note that the Windows check is additionally duplicated -- deliberately and
# with a pointer back here -- as an early guard in the top-level CMakeLists.txt
# ahead of project(). See the comment there for why.
#------------------------------------------------------------------------------

# --- MSVC ---------------------------------------------------------------------
# Several spellings of the same release, because each consumer needs a different one:
#   OPENMS_MIN_MSVC_VERSION   - the MSVC_VERSION integer, for the check below
#   OPENMS_MIN_MSVC_TOOLSET   - the toolset number, as reported by %VCToolsVersion%
#   OPENMS_MIN_MSVC_IDE       - the VS release, which is what users actually install
#   OPENMS_MIN_MSVC_YEAR      - the VS product year, as it appears in the product name
#                               and in the default install path
#   OPENMS_MIN_MSVC_GENERATOR - the matching CMake generator name
#
# The binding constraint is Apache Arrow. Toolsets older than
# 14.44 carry a bug in the standard library's <chrono> formatter for zoned_time
# durations coarser than seconds (LWG-4124, fixed by microsoft/STL#5155), which
# Arrow's std::chrono backend trips over while compiling strftime for Date32.
# The failure surfaces as an unreadable template error inside <chrono>, hence
# the explicit guard.
set(OPENMS_MIN_MSVC_VERSION 1944)
set(OPENMS_MIN_MSVC_TOOLSET "14.44")
set(OPENMS_MIN_MSVC_IDE "17.14")
set(OPENMS_MIN_MSVC_YEAR "2022")
# The generator name embeds the VS major version (the "17" of 17.14) and the product
# year, so it changes whenever either does. Derived rather than spelled out, so there
# is still only one place to edit.
string(REGEX REPLACE "\\..*$" "" _openms_msvc_major "${OPENMS_MIN_MSVC_IDE}")
set(OPENMS_MIN_MSVC_GENERATOR "Visual Studio ${_openms_msvc_major} ${OPENMS_MIN_MSVC_YEAR}")
unset(_openms_msvc_major)

# --- GCC / Clang / AppleClang -------------------------------------------------
# Driven by C++23 feature support.
set(OPENMS_MIN_GCC "13")
set(OPENMS_MIN_CLANG "17")
# AppleClang reports its own version scheme; 16 is what ships with Xcode 16.
set(OPENMS_MIN_APPLECLANG "16")
set(OPENMS_MIN_XCODE "16")

# --- CMake --------------------------------------------------------------------
# Not a compiler, but the same "documentation must not drift from the build
# system" problem, and the reason for the sync comment on cmake_minimum_required()
# in the top-level CMakeLists.txt. Read it back from CMake itself rather than
# repeating the number, so it can never disagree.
#
# This is captured at top-level scope, which is also the only scope that declares
# a cmake_minimum_required() -- see the note there on why subdirectories must not
# declare their own.
set(OPENMS_MIN_CMAKE "${CMAKE_MINIMUM_REQUIRED_VERSION}")

#------------------------------------------------------------------------------
# Enforcement
#------------------------------------------------------------------------------
# Test CMAKE_CXX_COMPILER_ID rather than the MSVC variable: MSVC is also true
# for clang-cl, which has its own version scheme and is covered by the Clang
# branch below.
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  if(MSVC_VERSION LESS OPENMS_MIN_MSVC_VERSION)
    message(FATAL_ERROR
      "OpenMS requires Visual Studio ${OPENMS_MIN_MSVC_YEAR} version ${OPENMS_MIN_MSVC_IDE} or later "
      "(MSVC toolset ${OPENMS_MIN_MSVC_TOOLSET}+, MSVC_VERSION >= ${OPENMS_MIN_MSVC_VERSION}), "
      "but found MSVC_VERSION ${MSVC_VERSION} (cl.exe ${CMAKE_CXX_COMPILER_VERSION}).\n"
      "\n"
      "Please update Visual Studio.")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS OPENMS_MIN_GCC)
    message(FATAL_ERROR
      "OpenMS requires GCC ${OPENMS_MIN_GCC} or later for C++23 support, "
      "but found ${CMAKE_CXX_COMPILER_VERSION}.")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS OPENMS_MIN_APPLECLANG)
    message(FATAL_ERROR
      "OpenMS requires Apple Clang ${OPENMS_MIN_APPLECLANG} or later "
      "(ships with Xcode ${OPENMS_MIN_XCODE}) for C++23 support, "
      "but found ${CMAKE_CXX_COMPILER_VERSION}.")
  endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS OPENMS_MIN_CLANG)
    message(FATAL_ERROR
      "OpenMS requires Clang ${OPENMS_MIN_CLANG} or later for C++23 support, "
      "but found ${CMAKE_CXX_COMPILER_VERSION}.")
  endif()
else()
  message(STATUS
    "Compiler '${CMAKE_CXX_COMPILER_ID}' (${CMAKE_CXX_COMPILER_VERSION}) is not one of the "
    "officially supported compilers. No minimum version is enforced; the build may or may not work.")
endif()

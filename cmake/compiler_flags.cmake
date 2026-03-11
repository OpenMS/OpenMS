# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------


#------------------------------------------------------------------------------
# This cmake file handles all the project specific compiler flags
# It defines variables and options and provides functions to apply flags to targets

#------------------------------------------------------------------------------
# PART 1: Define variables and options
#------------------------------------------------------------------------------

# Custom compile flags
if (MY_CXX_FLAGS) ## do not change this name! it's used in configh.cmake
  message(STATUS "Custom compile flags: '${MY_CXX_FLAGS}' will be added to targets")
endif()

# SIMD extensions
set(x64_CPU "x86|AMD64") ## CMake returns 'x86-64' on Linux and 'AMD64' on Windows..
message(STATUS "Processor is : ${CMAKE_SYSTEM_PROCESSOR}")

# ARM processors: SIMDe will do the right thing upon detecting ARM
# https://github.com/simd-everywhere/simde/blob/master/simde/simde-arch.h#L117
# (neon instructions compile without error even if no compile flag is given -- as opposed to x64 intrinsics)

# Compiler-specific options
if(CMAKE_COMPILER_IS_GNUCXX)
  # GCC options
  option(ENABLE_GCC_WERROR "Enable -Werror on gcc compilers" OFF)
  if(ENABLE_GCC_WERROR)
    message(STATUS "Enable -Werror for gcc - note that this may not work on all compilers and system settings!")
  endif()
elseif(MSVC)
  # MSVC options
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # Clang options
  set(CMAKE_COMPILER_IS_CLANG true CACHE INTERNAL "Is CLang compiler (clang++)")
else()
  # Intel compiler options
  set(CMAKE_COMPILER_IS_INTELCXX true CACHE INTERNAL "Is Intel C++ compiler (icpc)")
endif()

# Platform-dependent options
include(CheckCXXCompilerFlag)
if(NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
  CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
  if(WITH_FPIC)
    message(STATUS "Position-independent code (-fPIC) is supported and will be applied to targets")
  endif()
endif()

# Conversion warnings
set(CXX_WARN_CONVERSION OFF CACHE BOOL "Enables warnings for type conversion problems (GCC only)")
message(STATUS "Compiler checks for conversion: ${CXX_WARN_CONVERSION}")

#------------------------------------------------------------------------------
# PART 2: Functions to apply compiler flags to targets
#------------------------------------------------------------------------------

# Function to add compiler flags to a target with proper PUBLIC/PRIVATE visibility
function(openms_add_compiler_flags target_name)
  #------------------------------------------------------------------------------
  # PUBLIC flags (propagated to dependent targets)
  #------------------------------------------------------------------------------
  
  # Language standard
  target_compile_features(${target_name} PUBLIC cxx_std_20)
  
  # Position-independent code
  if(NOT WIN32 AND WITH_FPIC)
    target_compile_options(${target_name} PUBLIC -fPIC)
  endif()
  
  # Essential preprocessor definitions
  if(MSVC)
    target_compile_definitions(${target_name} PUBLIC
      NOMINMAX  # coinor windows.h include bug workaround
    )
  endif()
  
  # SIMD extensions (PUBLIC for binary compatibility)
  # MSVC x64 defaults to SSE2 (128-bit); do NOT use /arch:AVX here because
  # AVX's 256-bit reductions change Eigen's floating-point evaluation order
  # vs the 128-bit paths used by GCC/Clang, causing cross-platform numeric
  # divergence in sensitive algorithms (e.g. normalizedCrossCorrelation).
  if(NOT MSVC AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "${x64_CPU}")
    target_compile_options(${target_name} PUBLIC -mssse3)
  endif()
  
  #------------------------------------------------------------------------------
  # PRIVATE flags (not propagated to dependent targets)
  #------------------------------------------------------------------------------
  
  # Warning controls
  if(CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(${target_name} PRIVATE
      -Wall -Wextra
      -Wno-unknown-pragmas
      -ffp-contract=off
      -Wno-unused-function
      -Wno-psabi
    )
    
    if(ENABLE_GCC_WERROR)
      target_compile_options(${target_name} PRIVATE -Werror)
    endif()
    
    if(CMAKE_GENERATOR STREQUAL "Eclipse CDT4 - Unix Makefiles")
      target_compile_options(${target_name} PRIVATE -fmessage-length=0)
    endif()
  elseif(MSVC)
    target_compile_options(${target_name} PRIVATE
      /wd4251 /wd4275  # disable dll-interface warning
      /wd4996          # disable deprecated functions warning
      /wd4661          # disable explicit template instantiation request warning
      /wd4503          # disable decorated name length exceeded warning9
      /bigobj          # for large object files
      /MP              # use multiple CPU cores
    )
    
    target_compile_definitions(${target_name} PRIVATE
      _SCL_SECURE_NO_WARNINGS
      _CRT_SECURE_NO_WARNINGS
      _CRT_SECURE_NO_DEPRECATE
      OPENMS_XERCESDLL  # xerces bug workaround
    )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${target_name} PRIVATE
      -ffp-contract=off
      -Wall -Wextra
      -Wno-sign-conversion
      -Wno-global-constructors
      -Wno-exit-time-destructors
      -Wno-weak-vtables
      -Wno-documentation-unknown-command
      -Wno-undef
      -Wno-documentation
      -Wno-source-uses-openmp
      -Wno-old-style-cast
      -Wno-unknown-warning-option
      -Wno-double-promotion
      -Wno-unused-template
      -Wno-conversion
      -Wno-float-equal
      -Wno-switch-enum
      -Wno-missing-prototypes
      -Wno-missing-variable-declarations
      -Wno-deprecated
      -Wno-covered-switch-default
      -Wno-date-time
      -Wno-missing-noreturn
    )
  endif()
  
  # Conversion warnings
  if(CXX_WARN_CONVERSION AND CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(${target_name} PRIVATE -Wconversion)
  endif()
  
  # Custom flags
  if(MY_CXX_FLAGS)
    separate_arguments(MY_CXX_FLAGS_LIST UNIX_COMMAND "${MY_CXX_FLAGS}")
    target_compile_options(${target_name} PRIVATE ${MY_CXX_FLAGS_LIST})
  endif()
endfunction()

# Function to add library-specific compiler flags
function(openms_add_library_compiler_flags target_name)
  # Add common flags first (both PUBLIC and PRIVATE)
  openms_add_compiler_flags(${target_name})
  
  # Library-specific PRIVATE flags
  if(CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(${target_name} PRIVATE -Wno-non-virtual-dtor)
  endif()
  
  # Address sanitizer (PRIVATE)
  if(ADDRESS_SANITIZER)
    add_asan_to_target(${target_name})
  endif()
endfunction()

# Function to add executable-specific compiler flags
function(openms_add_executable_compiler_flags target_name)
  # Add common flags first (both PUBLIC and PRIVATE)
  openms_add_compiler_flags(${target_name})
  
  # Executable-specific PRIVATE flags can be added here
  
  # Address sanitizer (PRIVATE)
  if(ADDRESS_SANITIZER)
    add_asan_to_target(${target_name})
  endif()
endfunction()

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: OpenMS Team $
# $Authors: OpenMS Team $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# This cmake file handles all the target-specific compiler flags
# It replaces the global compiler flags previously set in compiler_flags.cmake

# Function to add common compiler flags to a target
function(openms_add_common_compiler_flags target_name)
  if(CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(${target_name} PRIVATE
      -Wall -Wextra
      -Wno-unknown-pragmas
      -Wno-long-long 
      -Wno-unused-function
      -Wno-variadic-macros
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
      /wd4503          # disable decorated name length exceeded warning
      /wd4068          # disable unknown pragma warning
      /bigobj          # for large object files
      /MP              # use multiple CPU cores
    )
    
    target_compile_definitions(${target_name} PRIVATE
      _SCL_SECURE_NO_WARNINGS
      _CRT_SECURE_NO_WARNINGS
      _CRT_SECURE_NO_DEPRECATE
      OPENMS_XERCESDLL  # xerces bug workaround
      NOMINMAX          # coinor windows.h include bug workaround
    )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${target_name} PRIVATE
      -ffp-contract=off
      -Wall -Wextra
      -Wno-sign-conversion
      -Wno-long-long
      -Wno-padded
      -Wno-global-constructors
      -Wno-exit-time-destructors
      -Wno-weak-vtables
      -Wno-documentation-unknown-command
      -Wno-undef
      -Wno-documentation
      -Wno-source-uses-openmp
      -Wno-old-style-cast
      -Wno-c++98-compat
      -Wno-c++98-compat-pedantic
      -Wno-unknown-warning-option
      -Wno-double-promotion
      -Wno-unused-template
      -Wno-conversion
      -Wno-float-equal
      -Wno-switch-enum
      -Wno-missing-prototypes
      -Wno-missing-variable-declarations
      -Wno-deprecated
      -Wno-deprecated-register
      -Wno-covered-switch-default
      -Wno-date-time
      -Wno-missing-noreturn
    )
  endif()
  
  # Platform-dependent flags
  if(NOT WIN32 AND WITH_FPIC)
    target_compile_options(${target_name} PRIVATE -fPIC)
  endif()
  
  # Conversion warnings
  if(CXX_WARN_CONVERSION AND CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(${target_name} PRIVATE -Wconversion)
  endif()
  
  # SIMD extensions
  if(MSVC AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "${x64_CPU}")
    target_compile_options(${target_name} PRIVATE /arch:AVX)
  elseif(NOT MSVC AND ${CMAKE_SYSTEM_PROCESSOR} MATCHES "${x64_CPU}")
    target_compile_options(${target_name} PRIVATE -mssse3)
  endif()
  
  # Custom flags
  if(MY_CXX_FLAGS)
    separate_arguments(MY_CXX_FLAGS_LIST UNIX_COMMAND "${MY_CXX_FLAGS}")
    target_compile_options(${target_name} PRIVATE ${MY_CXX_FLAGS_LIST})
  endif()
endfunction()

# Function to add library-specific compiler flags
function(openms_add_library_compiler_flags target_name)
  # Add common flags first
  openms_add_common_compiler_flags(${target_name})
  
  # Library-specific flags
  if(CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(${target_name} PRIVATE -Wno-non-virtual-dtor)
  endif()
  
  # Address sanitizer
  if(ADDRESS_SANITIZER)
    add_asan_to_target(${target_name})
  endif()
endfunction()

# Function to add executable-specific compiler flags
function(openms_add_executable_compiler_flags target_name)
  # Add common flags first
  openms_add_common_compiler_flags(${target_name})
  
  # Executable-specific flags can be added here
  
  # Address sanitizer
  if(ADDRESS_SANITIZER)
    add_asan_to_target(${target_name})
  endif()
endfunction()
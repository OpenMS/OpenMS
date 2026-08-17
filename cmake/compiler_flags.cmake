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
  target_compile_features(${target_name} PUBLIC cxx_std_23)
  
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
  
  # C++23 <stdfloat> workaround for ARM: std::float16_t/float32_t/float64_t
  # conflict with ARM NEON typedefs of the same name when 'using namespace std;'
  # is active. Suppress only the three macros that conflict with NEON typedefs.
  # Do NOT suppress __STDCPP_FLOAT128_T__ or __STDCPP_BFLOAT16_T__ as those
  # have no NEON conflict and are needed by Boost.Math.
  if(NOT MSVC AND CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
    target_compile_options(${target_name} PUBLIC
      -U__STDCPP_FLOAT16_T__ -U__STDCPP_FLOAT32_T__ -U__STDCPP_FLOAT64_T__
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

  # When full OpenMP isn't available (typically Apple clang without libomp),
  # multithreading.cmake sets OPENMS_OMP_SIMD_FALLBACK_FLAG to -fopenmp-simd
  # so #pragma omp simd vectorization still works inside OpenMS. Applied PRIVATE:
  # no installed public OpenMS header contains `#pragma omp simd`, so downstream
  # consumers don't need this flag. Propagating it PUBLIC could break downstream
  # builds whose compiler doesn't accept -fopenmp-simd.
  # Not gated on architecture: -fopenmp-simd is useful on ARM too (NEON).
  if(DEFINED OPENMS_OMP_SIMD_FALLBACK_FLAG)
    target_compile_options(${target_name} PRIVATE ${OPENMS_OMP_SIMD_FALLBACK_FLAG})
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

#------------------------------------------------------------------------------
# Keep symbols that come from static archives out of an executable's dynamic
# symbol table.
#
# Rationale: with ARROW_USE_STATIC=ON (the default) the Arrow/Parquet archives -
# including the Thrift code bundled inside them - are linked into libOpenMS.so
# *and* into every test executable that links Arrow/Parquet directly. The
# resulting weak template symbols (e.g. the static
# 'std::map<int, const char*>' enum-name tables built from
# apache::thrift::TEnumIterator) then exist in both binaries and are exported
# by both. Since the executable is searched first when the dynamic linker
# resolves global symbols, libOpenMS.so binds to the *executable's* copy and
# registers an atexit destructor for it. The object is consequently destroyed
# twice at process teardown - once via the executable's exit handlers and once
# via __cxa_finalize() when libOpenMS.so is unloaded - which corrupts the heap
# and aborts the process with "malloc_consolidate(): invalid chunk size" after
# all tests have already passed.
#
# Hiding the archive symbols makes each binary use its own copy, so every
# static object is destroyed exactly once. Executables never need to export
# symbols (nothing dlopen()s a test), so this is safe.
#
# GNU ld and lld only; the Apple linker has no equivalent option, and Windows
# does not have symbol interposition to begin with.
function(openms_hide_static_archive_symbols target_name)
  if(NOT WIN32 AND NOT APPLE)
    target_link_options(${target_name} PRIVATE "LINKER:--exclude-libs,ALL")
  endif()
endfunction()

#------------------------------------------------------------------------------
# Link feature 'needed': record a shared library as a dependency even if nothing
# has referenced it yet at the point where it appears on the link line.
#
# Most Linux distributions default their linker to --as-needed, which drops a
# shared library unless some object seen *earlier* already needs a symbol from
# it. That turns the order of target_link_libraries() arguments into a silent
# correctness requirement: a class test whose own object references nothing from
# libOpenMS keeps it only because libOpenMSTestFramework.a -- which does
# reference it -- happens to be listed first.
#
#   target_link_libraries(mytest OpenMSTestFramework "$<LINK_LIBRARY:needed,OpenMS>")
#
# states that requirement instead of encoding it in the argument order, so
# reordering the arguments can no longer break the link.
#
# The Apple linker resolves archives repeatedly and has no --as-needed, and
# Windows has no equivalent concept, so there the feature passes the library
# through unchanged. It is defined (not just supported) on every platform, so
# the generator expression stays valid everywhere.
if(NOT WIN32 AND NOT APPLE)
  set(CMAKE_CXX_LINK_LIBRARY_USING_needed
      "LINKER:--push-state,--no-as-needed"
      "<LINK_ITEM>"
      "LINKER:--pop-state")
else()
  set(CMAKE_CXX_LINK_LIBRARY_USING_needed "<LINK_ITEM>")
endif()
set(CMAKE_CXX_LINK_LIBRARY_USING_needed_SUPPORTED TRUE)

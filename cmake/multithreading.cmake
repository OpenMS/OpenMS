# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# OpenMP
#------------------------------------------------------------------------------
message(STATUS "OpenMP support requested: ${MT_ENABLE_OPENMP}")

if (MT_ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS CXX)

  if (NOT OPENMP_FOUND)
    # Full OpenMP runtime not found (typical on macOS with stock Apple clang).
    # Fall back to -fopenmp-simd, which enables `#pragma omp simd` vectorization
    # WITHOUT requiring libgomp/libomp at link time and WITHOUT defining _OPENMP.
    # Existing OpenMP::OpenMP_CXX consumers are gated by OPENMP_FOUND, so this
    # fallback does not activate any OpenMP-runtime-dependent code path.
    # Note: CHECK_CXX_COMPILER_FLAG caches OPENMS_HAS_FOPENMP_SIMD as an
    # INTERNAL cache variable (standard CMake behavior). If you upgrade the
    # compiler and want the probe to re-run, delete that cache entry or pass
    # -UOPENMS_HAS_FOPENMP_SIMD on the next configure.
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-fopenmp-simd" OPENMS_HAS_FOPENMP_SIMD)
    if (OPENMS_HAS_FOPENMP_SIMD)
      # Plain (non-CACHE) variable on purpose: we want it to reset to undefined
      # on every reconfigure, so that if the user later installs libomp and
      # find_package(OpenMP) starts succeeding, the fallback flag stops being
      # applied. A CACHE INTERNAL value would persist and silently combine with
      # full OpenMP on subsequent reconfigures.
      set(OPENMS_OMP_SIMD_FALLBACK_FLAG "-fopenmp-simd")
      message(STATUS
        "Full OpenMP not found; falling back to -fopenmp-simd. "
        "#pragma omp simd will vectorize, but #pragma omp parallel will not. "
        "On macOS, install Homebrew's libomp (`brew install libomp`) and "
        "reconfigure with -DOpenMP_ROOT=$(brew --prefix libomp) for full OpenMP. "
        "On Linux this fallback is uncommon (libgomp/libomp ship with GCC and "
        "Clang); if it fires here, check that your toolchain installation is "
        "complete and that find_package(OpenMP) can locate it.")
    else()
      message(FATAL_ERROR
        "OpenMP was requested (MT_ENABLE_OPENMP=ON) but neither full OpenMP "
        "nor the -fopenmp-simd fallback is supported by this compiler. "
        "Pass -DMT_ENABLE_OPENMP=OFF to skip OpenMP entirely.")
    endif()
  endif()
endif()

if (OPENMP_FOUND)
  set(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE) # will install the MSVC OpenMP runtime libraries
  
  # MSVC requires /openmp:experimental to support #pragma omp simd directives.
  # The /openmp:experimental flag is a superset of /openmp, so it replaces
  # rather than appends to the standard OpenMP flag.
  if (MSVC)
    set(OpenMP_CXX_FLAGS "/openmp:experimental")
    # Update the imported target to use the experimental flag instead of /openmp
    if (TARGET OpenMP::OpenMP_CXX)
      set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "/openmp:experimental")
    endif()
  endif()
endif()

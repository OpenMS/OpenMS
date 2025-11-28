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
    message(FATAL_ERROR "OpenMP was requested (MT_ENABLE_OPENMP=ON) but not found. "
                        "Please install OpenMP support or disable it with -DMT_ENABLE_OPENMP=OFF")
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

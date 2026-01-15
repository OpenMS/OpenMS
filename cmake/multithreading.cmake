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
endif()

if (OPENMP_FOUND)
  set(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE) # will install the MSVC OpenMP runtime libraries
endif()

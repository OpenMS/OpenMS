# Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# This cmake file enables the STL debug mode

if (CMAKE_COMPILER_IS_GNUCXX)
	if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
		# add compiler flag
  	add_compile_options(/D_GLIBCXX_DEBUG)
  	message(STATUS "STL debug mode: ${STL_DEBUG}")
  else()
    message(WARNING "STL debug mode is supported for OpenMS debug mode only")
  endif()
else()
  message(WARNING "STL debug mode is supported for compiler GCC only")
endif()

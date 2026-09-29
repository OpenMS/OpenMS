# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------


#------------------------------------------------------------------------------
# This cmake file only handles the customization of internal path variables
# where the build system expects to find external libraries. Note that the
# actual libraries are found in the CMake files of the individual componenents.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# OPENMS_CONTRIB_LIBS: the prefix of a contrib build (OpenMS/contrib), searched
# before the rest of CMAKE_PREFIX_PATH to avoid mismatches with installed system
# libraries. The contrib is retired in favour of vcpkg (#10327): the option is
# deprecated and goes away after OpenMS 3.6. vcpkg builds get their search paths
# from the vcpkg toolchain, all other builds from CMAKE_PREFIX_PATH.
if(OPENMS_CONTRIB_LIBS)
  if(OPENMS_USE_VCPKG)
    message(FATAL_ERROR
      "OPENMS_CONTRIB_LIBS cannot be combined with OPENMS_USE_VCPKG=ON, which takes the "
      "dependencies from vcpkg. Remove OPENMS_CONTRIB_LIBS, or set OPENMS_USE_VCPKG=OFF "
      "to build against the contrib.")
  endif()
  message(DEPRECATION
    "OPENMS_CONTRIB_LIBS is deprecated and will be removed after OpenMS 3.6, together "
    "with the contrib. Build the dependencies with vcpkg (cmake --preset <platform>-release, "
    "see the vcpkg install guide), or install them with a package manager and add their "
    "prefix to CMAKE_PREFIX_PATH.")
  list(INSERT CMAKE_PREFIX_PATH 0 ${OPENMS_CONTRIB_LIBS})
  list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
  list(REMOVE_ITEM CMAKE_PREFIX_PATH "") # Remove empty entries
endif()

#------------------------------------------------------------------------------
# Ensure Qt includes it's libs as SYSTEM
set(QT_INCLUDE_DIRS_NO_SYSTEM Off)
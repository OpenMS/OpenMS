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
# OPENMS_CONTRIB_LIBS pointed the build at a contrib build (OpenMS/contrib). The
# option went away after OpenMS 3.6, together with the contrib (#10327). Stop
# instead of ignoring it: a build that silently dropped the prefix would pick up
# other libraries than the ones the user asked for.
if(DEFINED OPENMS_CONTRIB_LIBS)
  message(FATAL_ERROR
    "OPENMS_CONTRIB_LIBS was removed after OpenMS 3.6, together with the contrib. "
    "Add the prefix of the dependencies to CMAKE_PREFIX_PATH, or let vcpkg build them "
    "(cmake --preset <platform>-release, see 'cmake --list-presets'). "
    "A build directory that still has the option cached needs 'cmake -U OPENMS_CONTRIB_LIBS'.")
endif()

#------------------------------------------------------------------------------
# Ensure Qt includes it's libs as SYSTEM
set(QT_INCLUDE_DIRS_NO_SYSTEM Off)
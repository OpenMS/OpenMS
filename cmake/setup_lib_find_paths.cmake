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
# Specify the path to the contrib-build. Will be added to CMAKE_PREFIX_PATH
# for searching and as first entry in the includes/libraries to avoid
# mismatches with installed system libraries
if(NOT OPENMS_CONTRIB_LIBS)
  message("Note: OPENMS_CONTRIB_LIBS not set. Unless you are certain that you have all contributing libraries in system paths, please specify an explicit path to the built contrib libraries via
-DOPENMS_CONTRIB_LIBS")
else()
  list(INSERT CMAKE_PREFIX_PATH 0 ${OPENMS_CONTRIB_LIBS})
  list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
  list(REMOVE_ITEM CMAKE_PREFIX_PATH "") # Remove empty entries
endif()

#------------------------------------------------------------------------------
# Ensure Qt includes it's libs as SYSTEM
set(QT_INCLUDE_DIRS_NO_SYSTEM Off)

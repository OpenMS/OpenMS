# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------


set(CPACK_RPM_PACKAGE_VENDOR "OpenMS developers <open-ms-general@lists.sourceforge.net>")
if((DEFINED ENV{CPACK_PACKAGE_FILE_NAME}) AND (NOT "$ENV{CPACK_PACKAGE_FILE_NAME}" STREQUAL ""))
  set(CPACK_PACKAGE_FILE_NAME "$ENV{CPACK_PACKAGE_FILE_NAME}")
else()
  if (OPENMS_64BIT_ARCHITECTURE)
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}-RedHat-Linux-x86_64")
  else()
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}-RedHat-Linux-x86")
  endif()
endif()

## Debug for now.
set(CPACK_RPM_PACKAGE_DEBUG ON)

## TODO also install headers? make a dev package configuration?
set(CPACK_COMPONENTS_ALL applications doc library share ${THIRDPARTY_COMPONENT_GROUP})

SET(CPACK_RPM_PACKAGE_LICENSE "BSD clause 3")
#SET(CPACK_RPM_PACKAGE_AUTOREQ ON)
#SET(CPACK_RPM_PACKAGE_AUTOPROV ON)
SET(CPACK_RPM_PACKAGE_AUTOREQPROV ON)
SET(CPACK_RPM_COMPRESSION_TYPE xz)
SET(CPACK_RPM_INSTALL_WITH_EXEC ON)
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "package for LC/MS data management and analysis")
SET(CPACK_PACKAGE_DESCRIPTION "
 OpenMS is a package for LC/MS data management and analysis. OpenMS
 offers an infrastructure for the development of mass
 spectrometry-related software and powerful 2D and 3D visualization
 solutions.
 .
 TOPP (the OpenMS PiPeline) is a pipeline for the analysis
 of HPLC/MS data. It consists of a set of numerous small applications
 that can be chained together to create analysis pipelines tailored
 for a specific problem."
 )

## Create own target because you cannot "depend" on the internal target 'package'
add_custom_target(dist
  COMMAND cpack -G ${CPACK_GENERATOR}
  COMMENT "Building ${CPACK_GENERATOR} package"
)

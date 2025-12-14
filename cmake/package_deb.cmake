# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------


set(CPACK_DEBIAN_PACKAGE_MAINTAINER "OpenMS Inc <info@openms.de>")
if((DEFINED ENV{CPACK_PACKAGE_FILE_NAME}) AND (NOT "$ENV{CPACK_PACKAGE_FILE_NAME}" STREQUAL ""))
  set(CPACK_PACKAGE_FILE_NAME "$ENV{CPACK_PACKAGE_FILE_NAME}")
else()
  # Use CMAKE_SYSTEM_PROCESSOR to detect actual architecture
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}-Debian-Linux-${CMAKE_SYSTEM_PROCESSOR}")
endif()

## CPack issues when building the package.
## https://bugs.launchpad.net/ubuntu/+source/cmake/+bug/972419
## https://ubuntuforums.org/showthread.php?t=2316865
## Workaround after packaging: https://cmake.org/pipermail/cmake/2012-May/050483.html
## Following needs CMake 3.7+. Just install from cmake.org
set(CPACK_DEBIAN_ARCHIVE_TYPE "gnutar")

## We usually do not want to ship things like stdlib or glibc. Could mess up a system slighlty, when installed system wide
#include(InstallRequiredSystemLibraries)

# Don't add RPATH
SET(CMAKE_SKIP_INSTALL_RPATH TRUE)

## Try autogeneration of dependencies:
## This may result in non-standard package names in the dependencies (e.g. when using Qt from a Thirdparty repo)
## It also will add system dependencies like a minimum glibc or gomp version (not necessarily bad)
##set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

## Debug for now. Not much output.
set(CPACK_DEBIAN_PACKAGE_DEBUG ON)

## TODO also install headers? make a dev package configuration?
set(CPACK_COMPONENTS_ALL applications doc library share ${THIRDPARTY_COMPONENT_GROUP})

## TODO we only need to put dependencies on shared libs. But this depends on what is found and what is statically linked on build machine.
## We should probably use a full system-shared-libs-only machine for building. Then the deps should look similar to below.
#set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxerces-c-dev (>= 3.1.1), libeigen3-dev, libboost-dev (>= 1.54.0), libboost-iostreams-dev (>= 1.54.0), libboost-date-time-dev (>= 1.54.0), libboost-math-dev (>= 1.54.0), libsvm-dev (>= 3.12), libglpk-dev (>= 4.52.1), zlib1g-dev (>= 1.2.7), libbz2-dev (>= 1.0.6), libqt4-dev (>= 4.8.2), libqt4-opengl-dev (>= 4.8.2), libqtwebkit-dev (>= 2.2.1), coinor-libcoinutils-dev (>= 2.6.4)")

## Autogeneration with SHLIBDEPS will add to this variable. For now we include most things statically and require the standard Qt package and libc6 only.
## (only available in Ubuntu >=17.10). For older Ubuntu, dependencies can be installed from a thirdparty repo.
## External libraries from Debian repositories should be listed here to avoid file conflicts
set(CPACK_DEBIAN_PACKAGE_DEPENDS 
  "libqt6svg6 (>= 6.2.2), libc6 (>= 2.28), libqt6widgets6t64 (>= 6.2.2) | libqt6widgets6 (>= 6.2.2), libqt6gui6t64 (>= 6.2.2) | libqt6gui6 (>= 6.2.2), libqt6core6t64 (>= 6.2.2) | libqt6core6 (>= 6.2.2), libyaml-cpp0.7 | libyaml-cpp0.8, libsqlite3-0 (>= 3.35.0), libsqlitecpp0 (>= 3.0.0)")

SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
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

## TODO make postinstall script that sets OPENMS_DATA_PATH

# For source packages add build dependencies. Not used and not tested.
#set(CPACK_DEBIAN_PACKAGE_BUILDS_DEPENDS "debhelper (>= 9), dpkg-dev (>= 1.16.1~), cmake (>= 2.6.3), imagemagick, doxygen (>= 1.8.1.2), graphviz, texlive-extra-utils, texlive-latex-extra, latex-xcolor, texlive-font-utils, ghostscript, texlive-fonts-recommended"

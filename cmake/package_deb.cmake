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
set(CPACK_DEBIAN_ARCHIVE_TYPE "gnutar")

## We usually do not want to ship things like stdlib or glibc. Could mess up a system slighlty, when installed system wide
#include(InstallRequiredSystemLibraries)

## Some libraries reach the staging tree with empty RUNPATH entries, which the loader
## reads as the current working directory: the release workflow links the binaries
## under PACKAGE_TYPE=none and packages them after reconfiguring (the script explains
## how that leaves them behind). The script removes them before the package is built.
list(APPEND CPACK_PRE_BUILD_SCRIPTS "${CMAKE_CURRENT_LIST_DIR}/cpack_clean_runpath.cmake")

## Derive the dependencies from the built binaries; a hand-written list goes stale and
## installs on systems the binaries cannot run on. Depends is what this derives and
## nothing else. A library the package ships itself (SQLite, in lib/) needs no Debian
## package, and a hand-written "<pkg>t64 | <pkg>" alternative would not help: Depends
## entries are ANDed, and dpkg-shlibdeps derives the plain <pkg>t64 name.
##
## This was on once before (#10202) and had to come off again (#10207), because the
## staging tree carried foreign-architecture binaries: ThermoRawFileParser's NuGet
## runtimes/<rid>/ tree ships a libMono.Unix.so for seven runtime identifiers
## (android-arm, android-arm64, android-x64, android-x86, linux-arm, linux-arm64,
## linux-x64), and on any one host five of the seven are foreign -- see the binary
## list in the failing run 35434550390, where they are the only such files.
## dpkg-shlibdeps answers each with "cannot find library libc.so.6 needed by
## ... (ELF format: ...; abi: ...)", an error --ignore-missing-info does not cover, so
## CPack aborted before writing the package. install_thirdparty_folder()
## (cmake/install_macros.cmake) now installs only the runtime identifiers this build
## targets, so nothing foreign reaches the staging tree. The private libraries the tool
## still cannot resolve (libOpenMS_CLI.so and its siblings) are same-architecture, and
## those come back as warnings rather than errors, so they do not stop it.
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

## Debug for now. Not much output.
set(CPACK_DEBIAN_PACKAGE_DEBUG ON)

## TODO also install headers? make a dev package configuration?
## The libraries come in layers (cmake/install_macros.cmake): library (core),
## library_cli (TOPP tool framework, needed by the TOPP tools) and library_gui.
## 'Applications' as install_tool() registers it: CPack installs a component by
## name and a mismatch would package none of the TOPP tools.
set(CPACK_COMPONENTS_ALL Applications doc library library_cli share ${THIRDPARTY_COMPONENT_GROUP})
if(WITH_GUI)
  list(APPEND CPACK_COMPONENTS_ALL library_gui)
endif()

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
#set(CPACK_DEBIAN_PACKAGE_BUILDS_DEPENDS "debhelper (>= 9), dpkg-dev (>= 1.16.1~), cmake (>= 2.6.3), imagemagick, doxygen (>= 1.8.1.2), graphviz, texlive-extra-utils, texlive-latex-extra, latex-xcolor, texlive-font-utils, texlive-fonts-recommended"

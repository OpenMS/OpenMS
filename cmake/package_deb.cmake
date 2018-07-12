set(CPACK_DEBIAN_PACKAGE_MAINTAINER "OpenMS developers <open-ms-general@lists.sourceforge.net>")
if (OPENMS_64BIT_ARCHITECTURE)
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}-Debian-Linux-x86_64")
else()
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}-Debian-Linux-x86")
endif()
set(CPACK_GENERATOR "DEB")

## CPack issues when building the package.
## https://bugs.launchpad.net/ubuntu/+source/cmake/+bug/972419
## https://ubuntuforums.org/showthread.php?t=2316865
## Workaround after packaging: https://cmake.org/pipermail/cmake/2012-May/050483.html
## Following needs CMake 3.7+. Just install from cmake.org
set(CPACK_DEBIAN_ARCHIVE_TYPE "gnutar")

## Not sure if we should ship them. Could mess up a system slighlty, when installed system wide?
#include(InstallRequiredSystemLibraries)

## Try autogeneration of dependencies
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
## Debug for now. Not much output.
set(CPACK_DEBIAN_PACKAGE_DEBUG ON)

## TODO also install headers? make a dev package configuration?
set(CPACK_COMPONENTS_ALL applications doc library share ${THIRDPARTY_COMPONENT_GROUP})

## TODO we only need to put dependencies on shared libs. But depends on what is found and what is statically linked on build machine.
## We should probably use a full system-shared-libs-only machine for building.
#set(CPACK_DEBIAN_PACKAGE_DEPENDS "libsqlite3-dev, libxerces-c-dev (>= 3.1.1), libeigen3-dev, libwildmagic-dev, libboost-dev (>= 1.54.0), libboost-iostreams-dev (>= 1.54.0), libboost-date-time-dev (>= 1.54.0), libboost-math-dev (>= 1.54.0), libsvm-dev (>= 3.12), libglpk-dev (>= 4.52.1), zlib1g-dev (>= 1.2.7), libbz2-dev (>= 1.0.6), libqt4-dev (>= 4.8.2), libqt4-opengl-dev (>= 4.8.2), libqtwebkit-dev (>= 2.2.1), coinor-libcoinutils-dev (>= 2.6.4)")

SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "package for LC/MS data management and analysis")
SET(CPACK_PACKAGE_DESCRIPTION "
 OpenMS is a package for LC/MS data management and analysis. OpenMS
 offers an infrastructure for the development of mass
 spectrometry-related software and powerful 2D and 3D visualization
 solutions.
 .
 TOPP (the OpenMS proteomic pipeline) is a pipeline for the analysis
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
#set(CPACK_DEBIAN_PACKAGE_BUILDS_DEPENDS "debhelper (>= 9), dpkg-dev (>= 1.16.1~), cmake (>= 2.6.3), imagemagick, doxygen (>= 1.8.1.2), graphviz, seqan-dev (>= 1.4.1), texlive-extra-utils, texlive-latex-extra, latex-xcolor, texlive-font-utils, ghostscript, texlive-fonts-recommended"

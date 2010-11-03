set(INSTALL_FORCE_DOC TRUE) ## force documentation to be present
include(cmake/install_common.cmake)

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "bielow@mi.fu-berlin.de")
set(CPACK_GENERATOR "DEB")
set(CPACK_COMPONENTS_ALL applications library share)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "cmake (>= 2.6), g++ (>= 4.0), autoconf (>= 2.6), qt4-dev-tools, patch, libtool")

## should be the last include
include(CPack)

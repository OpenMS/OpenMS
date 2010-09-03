set(INSTALL_FORCE_DOC TRUE) ## force documentation to be present
include(cmake/install_common.cmake)

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "bielow@mi.fu-berlin.de")
set(CPACK_GENERATOR "DEB")
set(CPACK_COMPONENTS_ALL applications library share)

## should be the last include
include(CPack)

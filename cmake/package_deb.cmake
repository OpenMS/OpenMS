set(INSTALL_FORCE_DOC TRUE) ## force documentation to be present
include(cmake/install_common.cmake)

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "bielow@mi.fu-berlin.de")
set(CPACK_GENERATOR "DEB")
set(CPACK_COMPONENTS_ALL applications library share)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libqt4-core|libqt4core, libqt4-gui|libqt4gui, libqt4-network, libqt4-opengl, libqt4-sql, libqt4-svg, libqt4-webkit, libqt4-xmlpatterns")

## should be the last include
include(CPack)

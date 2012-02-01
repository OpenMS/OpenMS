set(INSTALL_FORCE_DOC TRUE) ## force documentation to be present
include(cmake/install_common.cmake)

if (OPENMS_64BIT_ARCHITECTURE)
	set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
endif()
set(CPACK_GENERATOR "RPM")
set(CPACK_RPM_PACKAGE_REQUIRES "qt")
set(CPACK_RPM_POST_INSTALL_SCRIPT_FILE "${CMAKE_SOURCE_DIR}/cmake/run_ldconfig.sh")
set(CPACK_COMPONENTS_ALL applications library share)

## should be the last include
include(CPack)

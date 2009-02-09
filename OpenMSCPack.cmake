set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_NAME "OpenMS")
set(CPACK_PACKAGE_VENDOR "OpenMS.de")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMS - A framework for mass spectrometry")
set(CPACK_PACKAGE_VERSION "1.3.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "3")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMS 1.3")

# TODO: do we need to double this stuff
install(TARGETS OpenMS 
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT libraries)

foreach(i ${UTILS_executables})
	install(TARGETS ${i} RUNTIME DESTINATION bin COMPONENT applications)
endforeach()

set(CPACK_COMPONENTS_ALL applications libraries headers share contrib)

# should be the last include
include(CPack)


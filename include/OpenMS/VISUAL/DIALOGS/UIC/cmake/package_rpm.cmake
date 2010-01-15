## variables we need to distinguish 32Bit/64Bit versions
set(OPENMS_LIB_INSTALL_PATH "lib")
#set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86")
if (OPENMS_64BIT_ARCHITECTURE)
	set(OPENMS_LIB_INSTALL_PATH "lib64")
	set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
endif()

## CPack installation and packaging procedures
install(TARGETS OpenMS 
  LIBRARY DESTINATION ${OPENMS_LIB_INSTALL_PATH}
  ARCHIVE DESTINATION ${OPENMS_LIB_INSTALL_PATH}
  COMPONENT library)

## install utils
foreach(util ${UTILS_executables})
  install(TARGETS ${util}
    RUNTIME DESTINATION bin
    BUNDLE DESTINATION bin
    COMPONENT applications)
endforeach()

## install TOPP Tools
foreach(TOPP_exe ${TOPP_executables})
  INSTALL(TARGETS ${TOPP_exe} 
    RUNTIME DESTINATION bin
    BUNDLE DESTINATION bin
    COMPONENT applications)
endforeach()

## install share
INSTALL(DIRECTORY share/
  DESTINATION share
  COMPONENT share
  REGEX ".svn" EXCLUDE)

## install the documentation and the tutorials
install(FILES     doc/index.html      DESTINATION share/OpenMS/doc COMPONENT doc)
install(DIRECTORY doc/html            DESTINATION share/OpenMS/doc COMPONENT doc REGEX ".svn" EXCLUDE)
install(FILES doc/OpenMS_tutorial.pdf DESTINATION share/OpenMS/doc COMPONENT doc)
install(FILES doc/TOPP_tutorial.pdf   DESTINATION share/OpenMS/doc COMPONENT doc)

set(CPACK_GENERATOR "RPM")
set(CPACK_RPM_PACKAGE_REQUIRES "Qt4")
set(CPACK_COMPONENTS_ALL applications library share)

## should be the last include
include(CPack)
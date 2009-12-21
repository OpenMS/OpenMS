if (WIN32)
	set(OPENMS_LIB_INSTALL_PATH "bin")
else()	## Linux & MacOS
	## variables we need to distinguish 32Bit/64Bit versions
	set(OPENMS_LIB_INSTALL_PATH "lib")
	if (OPENMS_64BIT_ARCHITECTURE)
		set(OPENMS_LIB_INSTALL_PATH "lib64")
	endif()
endif()
	
## CPack installation and packaging procedures
install(TARGETS OpenMS EXPORT OpenMSLibExportGroup
  LIBRARY DESTINATION ${OPENMS_LIB_INSTALL_PATH}
  ARCHIVE DESTINATION ${OPENMS_LIB_INSTALL_PATH}
	RUNTIME DESTINATION ${OPENMS_LIB_INSTALL_PATH}
  COMPONENT library)

## create script that allows external projects to use our OpenMS lib
install(EXPORT OpenMSLibExportGroup DESTINATION cmake/ FILE ${CF_LibOpenMSExport})
		  
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
INSTALL(DIRECTORY share/		# warning: that slash(/) is important here, otherwise the whole directory (not its content) will be copied!
  DESTINATION share
  COMPONENT share
  PATTERN ".svn" EXCLUDE)

## install the documentation and the tutorials
install(FILES     ${PROJECT_BINARY_DIR}/doc/index.html      DESTINATION share/OpenMS/doc COMPONENT doc)
install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html            DESTINATION share/OpenMS/doc COMPONENT doc PATTERN ".svn" EXCLUDE)
install(FILES ${PROJECT_BINARY_DIR}/doc/OpenMS_tutorial.pdf DESTINATION share/OpenMS/doc COMPONENT doc)
install(FILES ${PROJECT_BINARY_DIR}/doc/TOPP_tutorial.pdf   DESTINATION share/OpenMS/doc COMPONENT doc)


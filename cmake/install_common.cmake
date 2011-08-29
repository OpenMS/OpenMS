if (WIN32)
	set(OPENMS_LIB_INSTALL_PATH "bin")
else()	## Linux & MacOS
  set(OPENMS_LIB_INSTALL_PATH "lib")
endif()

## CPack installation and packaging procedures
install(TARGETS OpenMS EXPORT OpenMSLibExportGroup
  LIBRARY DESTINATION ${OPENMS_LIB_INSTALL_PATH}
  ARCHIVE DESTINATION ${OPENMS_LIB_INSTALL_PATH}
	RUNTIME DESTINATION ${OPENMS_LIB_INSTALL_PATH}
  COMPONENT library)

		  
## install UTILS tools (using 'install(TARGETS UTILS...)' as shortcut does not work)
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

## install GUI Tools
foreach(TOPP_exe ${GUI_executables})
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

  
if (INSTALL_FORCE_DOC OR DOXYGEN_FOUND)
  # 'doc' target exists
  install(CODE "MESSAGE(\"Installing documentation created from 'doc' target. Make sure this target was called before!\")")
  ## this does not work yet (CMake 'bug') :  add_dependencies(install doc) ## force 'doc' to be build before installing the resulting files
  install(FILES     ${PROJECT_BINARY_DIR}/doc/index.html      DESTINATION share/OpenMS/doc COMPONENT doc) 
  install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html            DESTINATION share/OpenMS/doc COMPONENT doc PATTERN ".svn" EXCLUDE) 
  if (INSTALL_FORCE_DOC OR DOC_TUTORIALS_ACTIVE)
    # 'doc_tutorials' target exists
    install(CODE "MESSAGE(\"Installing documentation created from 'doc_tutorials' target. Make sure this target was called before!\")")
    ## this does not work yet (CMake 'bug') :  add_dependencies(install doc_tutorials) ## force 'doc_tutorials' to be build before installing the resulting files
    install(FILES ${PROJECT_BINARY_DIR}/doc/OpenMS_tutorial.pdf DESTINATION share/OpenMS/doc COMPONENT doc) 
    install(FILES ${PROJECT_BINARY_DIR}/doc/TOPP_tutorial.pdf   DESTINATION share/OpenMS/doc COMPONENT doc) 
  else() 
    Message(STATUS "Latex missing. Disabling 'doc_tutorials' target installation!") 
  endif() 
else() 
  Message(STATUS "Doxygen missing. Disabling all documentation targets installation!") 
endif() 

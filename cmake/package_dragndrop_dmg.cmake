set(CPACK_GENERATOR "DragNDrop")

# the OpenMS library
install (TARGETS OpenMS
	LIBRARY DESTINATION OpenMS-1.4/lib
	ARCHIVE DESTINATION OpenMS-1.4/lib
	COMPONENT library)

set(CMAKE_INSTALL_RPATH "../../../lib")
# the TOPP tools
foreach(TOPP_exe ${TOPP_executables})
	install(TARGETS ${TOPP_exe}
		RUNTIME DESTINATION OpenMS-1.4/bin
		BUNDLE DESTINATION OpenMS-1.4/Applications/
		COMPONENT applications)
endforeach()

# share dir
install(DIRECTORY share/
	DESTINATION OpenMS-1.4/share
	COMPONENT share
	REGEX ".svn" EXCLUDE)

# install the documentation and the tutorials
install(FILES     doc/index.html      		DESTINATION OpenMS-1.4/ COMPONENT doc)
install(DIRECTORY doc/html            		DESTINATION OpenMS-1.4/ COMPONENT doc REGEX ".svn" EXCLUDE)
install(FILES 		doc/OpenMS_tutorial.pdf DESTINATION OpenMS-1.4/ COMPONENT doc)
install(FILES 		doc/TOPP_tutorial.pdf   DESTINATION OpenMS-1.4/ COMPONENT doc)

# install the TOPP command shell
install(FILES cmake/MacOSX/TOPP-shell.command	DESTINATION OpenMS-1.4/ PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ COMPONENT TOPPShell)
install(FILES cmake/MacOSX/TOPP_bash_profile  DESTINATION OpenMS-1.4/ RENAME .TOPP_bash_profile COMPONENT TOPPShell)

include(CPack)

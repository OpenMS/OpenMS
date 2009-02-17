set(CPACK_GENERATOR "PackageMaker")

# the OpenMS library
install (TARGETS OpenMS
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	COMPONENT library)

# the TOPP tools
foreach(TOPP_exe ${TOPP_executables})
	install(TARGETS ${TOPP_exe}
		RUNTIME DESTINATION bin
		BUNDLE DESTINATION /Applications/
		COMPONENT applications)
endforeach()

# share dir
install(DIRECTORY share/
	DESTINATION share
	COMPONENT share
	REGEX ".svn" EXCLUDE)

# install the documentation and the tutorials
install(FILES     doc/index.html      DESTINATION share/OpenMS/doc COMPONENT doc)
install(DIRECTORY doc/html            DESTINATION share/OpenMS/doc COMPONENT doc REGEX ".svn" EXCLUDE)
install(FILES doc/OpenMS_tutorial.pdf DESTINATION share/OpenMS/doc COMPONENT doc)
install(FILES doc/TOPP_tutorial.pdf   DESTINATION share/OpenMS/doc COMPONENT doc)


include(CPack)
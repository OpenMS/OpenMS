set(CPACK_GENERATOR "DragNDrop")

# this option to take effect requires at least version 2.6.4 of cmake, which is not released yet
# therefore only the developer version can be used to make background images
# The file should be accompanied with a file named background png, stored in the sam folder
set(CPACK_DMG_DS_STORE ${PROJECT_SOURCE_DIR}/cmake/MacOSX/DS_store)

# add start-up script for TOPPView.app into Bundle, and move real TOPPView executable to TOPPView.exe
# call scripts to fix issues with wrongly referenced Qt libs 
add_custom_command(TARGET TOPPView POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/MacOS/TOPPView ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/MacOS/TOPPView.exe 
		COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/MacOS/TOPPView 
		COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPPView ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/MacOS
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/install_qt.bash ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/lib ${CMAKE_INSTALL_NAME_TOOL}
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_bundles.bash ${CMAKE_INSTALL_NAME_TOOL} ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${PROJECT_BINARY_DIR}/lib
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_openms.bash ${CMAKE_INSTALL_NAME_TOOL} ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/lib
    )

# add start-up script for INIFileEditor.app into Bundle, and move real INIFileEditor executable to INIFileEditor.exe
add_custom_command(TARGET INIFileEditor POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/MacOS/INIFileEditor ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/MacOS/INIFileEditor.exe 
    COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/MacOS/INIFileEditor 
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/MacOSX/INIFileEditor ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/MacOS)

# the OpenMS library
install (TARGETS OpenMS
	LIBRARY DESTINATION OpenMS-1.4/lib
	ARCHIVE DESTINATION OpenMS-1.4/lib
	COMPONENT library)

# Qt libs, hack! Copy files to destination first
install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/ DESTINATION OpenMS-1.4/lib/ COMPONENT library REGEX "libOpenMS.[A-Za-z]*" EXCLUDE)

# the TOPP tools
foreach(TOPP_exe ${TOPP_executables})
  if(NOT ${TOPP_exe} STREQUAL "TOPPView" AND NOT ${TOPP_exe} STREQUAL "INIFileEditor")
    # call scripts to fix issues with wrongly referenced Qt libs 
    add_custom_command(TARGET ${TOPP_exe} POST_BUILD
      COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_TOPP.bash ${CMAKE_INSTALL_NAME_TOOL} ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TOPP_exe}
      )
  endif()
	
  install(TARGETS ${TOPP_exe}
		RUNTIME DESTINATION OpenMS-1.4/TOPP
		BUNDLE DESTINATION OpenMS-1.4/
		COMPONENT applications)
endforeach()

# share dir
install(DIRECTORY share/
	DESTINATION OpenMS-1.4/share
	COMPONENT share
	REGEX ".svn" EXCLUDE)

# install the documentation and the tutorials
install(FILES     ${PROJECT_BINARY_DIR}/doc/index.html      		DESTINATION OpenMS-1.4/Documentation/ RENAME OpenMSAndTOPPDocumentation.html COMPONENT doc)
install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html            		DESTINATION OpenMS-1.4/Documentation/ COMPONENT doc REGEX ".svn" EXCLUDE)
install(FILES 		${PROJECT_BINARY_DIR}/doc/OpenMS_tutorial.pdf DESTINATION OpenMS-1.4/Documentation/ COMPONENT doc)
install(FILES 		${PROJECT_BINARY_DIR}/doc/TOPP_tutorial.pdf   DESTINATION OpenMS-1.4/Documentation/ COMPONENT doc)

# install the TOPP command shell
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP-shell.command DESTINATION OpenMS-1.4/ RENAME TOPP-shell PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ COMPONENT TOPPShell)
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP_bash_profile  DESTINATION OpenMS-1.4/ RENAME .TOPP_bash_profile COMPONENT TOPPShell)

# install background image for the folder
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/background.png DESTINATION OpenMS-1.4/share/OpenMS/ COMPONENT share)

include(CPack)

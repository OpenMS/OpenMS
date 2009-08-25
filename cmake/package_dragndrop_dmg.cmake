set(CPACK_GENERATOR "DragNDrop")

# this option to take effect requires at least version 2.6.4 of cmake
# The file should be accompanied with a file named background png, stored in the same folder
# To create the file convert a sample dmg into a read/write one, select background image
# size ... and copy the .DS_store the file in the MacOSX dir
# does not work currently, file is simply copied (see below)
#set(CPACK_DMG_DS_STORE ${PROJECT_SOURCE_DIR}/cmake/MacOSX/DS_store)

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

# add start-up script for TOPPAS.app into Bundle, and move real TOPPAS executable to TOPPAS.exe
add_custom_command(TARGET TOPPAS POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/MacOS/TOPPAS ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/MacOS/TOPPAS.exe 
    COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/MacOS/TOPPAS 
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPPAS ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/MacOS)

# the OpenMS library
install (TARGETS OpenMS
	LIBRARY DESTINATION OpenMS-1.5/lib
	ARCHIVE DESTINATION OpenMS-1.5/lib
	COMPONENT library)

# Qt libs, hack! Copy files to destination first
install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/ DESTINATION OpenMS-1.5/lib/ COMPONENT library REGEX "libOpenMS.[A-Za-z]*" EXCLUDE)

# the TOPP tools
foreach(TOPP_exe ${TOPP_executables})
  if(NOT ${TOPP_exe} STREQUAL "TOPPView" AND NOT ${TOPP_exe} STREQUAL "INIFileEditor" AND NOT ${TOPP_exe} STREQUAL "TOPPAS")
    # call scripts to fix issues with wrongly referenced Qt libs 
    add_custom_command(TARGET ${TOPP_exe} POST_BUILD
      COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_TOPP.bash ${CMAKE_INSTALL_NAME_TOOL} ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TOPP_exe}
      )
  endif()
	
  install(TARGETS ${TOPP_exe}
		RUNTIME DESTINATION OpenMS-1.5/TOPP
		BUNDLE DESTINATION OpenMS-1.5/
		COMPONENT applications)
endforeach()

# share dir
install(DIRECTORY share/
	DESTINATION OpenMS-1.5/share
	COMPONENT share
	REGEX ".svn" EXCLUDE)

# install the documentation and the tutorials
install(FILES     ${PROJECT_BINARY_DIR}/doc/index.html      		DESTINATION OpenMS-1.5/Documentation/ RENAME OpenMSAndTOPPDocumentation.html COMPONENT doc)
install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html            		DESTINATION OpenMS-1.5/Documentation/ COMPONENT doc REGEX ".svn" EXCLUDE)
install(FILES 		${PROJECT_BINARY_DIR}/doc/OpenMS_tutorial.pdf DESTINATION OpenMS-1.5/Documentation/ COMPONENT doc)
install(FILES 		${PROJECT_BINARY_DIR}/doc/TOPP_tutorial.pdf   DESTINATION OpenMS-1.5/Documentation/ COMPONENT doc)

# install the TOPP command shell
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP-shell.command DESTINATION OpenMS-1.5/ RENAME TOPP-shell PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ COMPONENT TOPPShell)
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP_bash_profile  DESTINATION OpenMS-1.5/ RENAME .TOPP_bash_profile COMPONENT TOPPShell)

# install background image for the folder
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/background.png DESTINATION OpenMS-1.5/share/OpenMS/ COMPONENT share)

# copy DS_store file into folder
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/DS_store DESTINATION . RENAME .DS_store COMPONENT share)

include(CPack)

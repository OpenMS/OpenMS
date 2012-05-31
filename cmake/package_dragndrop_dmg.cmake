set(CPACK_GENERATOR "DragNDrop")

# this option to take effect requires at least version 2.6.4 of cmake
# The file should be accompanied with a file named background png, stored in the same folder
# To create the file convert a sample dmg into a read/write one, select background image
# size ... and copy the .DS_store the file in the MacOSX dir
# does not work currently, file is simply copied (see below)
#set(CPACK_DMG_DS_STORE ${PROJECT_SOURCE_DIR}/cmake/MacOSX/DS_store)

# call scripts to fix issues with wrongly referenced Qt libs 

# create qt conf
file(WRITE "${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/Resources/qt.conf" "")
file(WRITE "${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/Resources/qt.conf" "")
file(WRITE "${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/Resources/qt.conf" "")


# qtmenu.nib
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/Resources/qt_menu.nib")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${QT_LIBRARY_DIR}/QtGui.framework/Resources/qt_menu.nib" "${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/Resources/qt_menu.nib")

file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/Resources/qt_menu.nib")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${QT_LIBRARY_DIR}/QtGui.framework/Resources/qt_menu.nib" "${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/Resources/qt_menu.nib")

file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/Resources/qt_menu.nib")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${QT_LIBRARY_DIR}/QtGui.framework/Resources/qt_menu.nib" "${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/Resources/qt_menu.nib")


# add openms fix script also here
INSTALL(CODE "
  execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/MacOS/TOPPView )
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPPView-resources/TOPPView 
    ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/MacOS)

  execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/MacOS/INIFileEditor) 
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/MacOSX/INIFileEditor-resources/INIFileEditor
    ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/MacOS)

  execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/MacOS/TOPPAS) 
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPPAS-resources/TOPPAS
    ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/MacOS)
"
  COMPONENT applications
)

install(CODE "
	execute_process(COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/install_qt.bash ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/lib ${CMAKE_INSTALL_NAME_TOOL})
	execute_process(COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_openms.bash ${CMAKE_INSTALL_NAME_TOOL} ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/lib)
"
  COMPONENT Library-Fixes
)


#pseudo-app install hack
install(DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/ 
  COMPONENT applications
  PATTERN ".svn" EXCLUDE
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE  
)

# hack to avoid permission problem
install(PROGRAMS  ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPAS.app/Contents/MacOS/TOPPAS 
  COMPONENT applications
  DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/TOPPAS.app/Contents/MacOS/
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

install(DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/ 
  COMPONENT applications
	PATTERN ".svn" EXCLUDE
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE  
)
install(PROGRAMS  ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/TOPPView.app/Contents/MacOS/TOPPView 
  COMPONENT applications
  DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/TOPPView.app/Contents/MacOS/
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

install(DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/ 
  COMPONENT applications
 	PATTERN ".svn" EXCLUDE
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE  
)
install(PROGRAMS  ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/INIFileEditor.app/Contents/MacOS/INIFileEditor 
  COMPONENT applications
  DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/INIFileEditor.app/Contents/MacOS/
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

# the OpenMS library
install (TARGETS OpenMS OpenMS_GUI
	LIBRARY DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/lib
	ARCHIVE DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/lib
	COMPONENT library)

# Qt libs, hack! Copy files to destination first
install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/ 
  DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/lib/ 
  COMPONENT library 
  REGEX "libOpenMS.[A-Za-z]*" EXCLUDE)

# the TOPP tools
foreach(TOPP_exe ${TOPP_executables})
  # call scripts to fix issues with wrongly referenced Qt libs 
  add_custom_command(TARGET ${TOPP_exe} POST_BUILD
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_TOPP.bash ${CMAKE_INSTALL_NAME_TOOL} ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TOPP_exe}
    )
	
  install(TARGETS ${TOPP_exe}
		RUNTIME DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/TOPP
		BUNDLE DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/
		COMPONENT applications)
endforeach()

# the UTILS
foreach(UTIL_exe ${UTILS_executables})
	# call scripts to fix issues with wrongly referenced Qt libs 
  add_custom_command(TARGET ${UTIL_exe} POST_BUILD
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_TOPP.bash ${CMAKE_INSTALL_NAME_TOOL} ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${UTIL_exe}
    )
  
	install(TARGETS ${UTIL_exe}
    RUNTIME DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/TOPP
    BUNDLE DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/
    COMPONENT applications)	
endforeach(UTIL_exe ${UTILS_executables})

# the GUI tools
foreach(GUI_exe ${GUI_executables})
  # call scripts to fix issues with wrongly referenced Qt libs 
  add_custom_command(TARGET ${GUI_exe} POST_BUILD
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_TOPP.bash ${CMAKE_INSTALL_NAME_TOOL} ${QT_LIBRARY_DIR} ${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${GUI_exe}
    )
	
  install(TARGETS ${GUI_exe}
		RUNTIME DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/TOPP
		BUNDLE DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/
		COMPONENT applications)
endforeach(GUI_exe ${GUI_executables})

# share dir
install(DIRECTORY share/
	DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/share
	COMPONENT share
	PATTERN ".svn" EXCLUDE)

# install the documentation and the tutorials
install(FILES     ${PROJECT_BINARY_DIR}/doc/index.html      		DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/ RENAME OpenMSAndTOPPDocumentation.html COMPONENT doc)
install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html            		DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/ COMPONENT doc PATTERN ".svn" EXCLUDE)
install(FILES 		${PROJECT_BINARY_DIR}/doc/OpenMS_tutorial.pdf DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/ COMPONENT doc)
install(FILES 		${PROJECT_BINARY_DIR}/doc/TOPP_tutorial.pdf   DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/ COMPONENT doc)

# install the TOPP command shell
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP-shell.command DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/ RENAME TOPP-shell PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ COMPONENT TOPPShell)
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP_bash_profile  DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/ RENAME .TOPP_bash_profile COMPONENT TOPPShell)

# install background image for the folder
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/background.png DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/share/OpenMS/ COMPONENT share)

# we fix this later by an external script
# copy DS_store file into folder
# install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/DS_store DESTINATION . RENAME .DS_store COMPONENT share)

include(CPack)

## fix the background image of the generated dmg
add_custom_target(final_package
  COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fixdmg.sh
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  COMMENT "Finalizing dmg image"
  DEPENDS dmg)
  
add_custom_target(dmg
  COMMAND cpack -G DragNDrop
  COMMENT "Building intermediate dmg package"
)

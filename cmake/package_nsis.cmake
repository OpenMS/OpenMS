## Windows installer

set(PACKAGE_LIB_DIR bin)
set(PACKAGE_BIN_DIR TOPP)

if (CMAKE_GENERATOR MATCHES ".*Win64.*")
  set(PLATFORM "64")
  set(ARCH "x64")
else()
  set(PLATFORM "32")
  set(ARCH "x86")
endif()

## Find redistributable to be installed by NSIS
if (NOT VC_REDIST_PATH)
	if(CMAKE_GENERATOR MATCHES ".*Visual Studio 1[1-9].*")
	  set(VC_REDIST_PATH "$ENV{VSINSTALLDIR}VC\\redist\\1033\\vcredist_${ARCH}.exe")
	else()
	  message(FATAL_ERROR "Variable VC_REDIST_PATH missing. Before Visual Studio 2012 you have to provide the file and its path on your own.")
	endif()
endif()
##TODO use following instead
# ########################################################### System runtime libraries
# set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP TRUE)
# include(InstallRequiredSystemLibraries)
# install(PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
#         DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/${PACKAGE_LIB_DIR}/
#         COMPONENT library)

configure_file(${PROJECT_SOURCE_DIR}/cmake/Windows/Cfg_Settings.nsh.in ${PROJECT_BINARY_DIR}/Cfg_Settings.nsh)

set(CPACK_GENERATOR NSIS)

set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}_Win${PLATFORM}")
set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS.ico")

## For now we fully  rely only on our NSIS template. Later we could use CMake install commands
## and the following to let CMake generate snippets for the NSIS script

# set(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS.ico")
# set(CPACK_NSIS_MUI_UNIICON "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS.ico")
# set(CPACK_NSIS_HELP_LINK "https://www.openms.de/getting-started")
# set(CPACK_NSIS_URL_INFO_ABOUT "https://www.openms.de")
# set(CPACK_NSIS_CONTACT "open-ms-general@lists.sourceforge.net")
# set(CPACK_NSIS_MENU_LINKS
#     "https://www.openms.de" "OpenMS Web Site")




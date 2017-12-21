## Windows installer
## TODO readd UAC plugin.. otherwise fails as normal user "Cant write file"
## No description of the search engines (subsections of Thirdparty)

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
	  set(VC_REDIST_PATH "$ENV{VSINSTALLDIR}VC\\redist\\1033")
	  set(VC_REDIST_EXE "vcredist_${ARCH}.exe")
	else()
	  message(FATAL_ERROR "Variable VC_REDIST_PATH missing."
	  "Before Visual Studio 2012 you have to provide the path"
	  "to the redistributable package of the VS you are using on your own.")
	endif()
endif()

##TODO try following instead once CMake generates NSIS commands for us.
# ########################################################### System runtime libraries
# set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP TRUE)
# include(InstallRequiredSystemLibraries)
# install(PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
#         DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/${PACKAGE_LIB_DIR}/
#         COMPONENT library)

configure_file(${PROJECT_SOURCE_DIR}/cmake/Windows/Cfg_Settings.nsh.in ${PROJECT_BINARY_DIR}/Cfg_Settings.nsh)

set(CPACK_GENERATOR NSIS)

set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-Win${PLATFORM}")
set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS.ico")

## Create own target because you cannot "depend" on the internal target 'package'
add_custom_target(dist
  COMMAND cpack -G ${CPACK_GENERATOR}
  COMMENT "Building ${CPACK_GENERATOR} package"
)

## TODO maybe find signtool and maybe check existence of ID in the beginning.
## ID needs to be installed beforehand. Rightclick a p12 file that has a cert for codesigning.
if (DEFINED SIGNING_IDENTITY AND NOT "${SIGNING_IDENTITY}" STREQUAL "") 
	add_custom_target(signed_dist
	  COMMAND signtool sign /v /n "${SIGNING_IDENTITY}" /t http://timestamp.digicert.com ${CPACK_PACKAGE_FILE_NAME}.exe
	  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	  COMMENT "Signing ${CPACK_PACKAGE_FILE_NAME}.exe with '${SIGNING_IDENTITY}'"
	  DEPENDS dist
	)
endif()

## For now we fully rely only on our NSIS template. Later we could use
## the following to let CMake generate snippets for the NSIS script
## Plus an additional entry in the nsis template (see CPack-NSIS docu)

# set(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS.ico")
# set(CPACK_NSIS_MUI_UNIICON "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS.ico")
# set(CPACK_NSIS_HELP_LINK "https://www.openms.de/getting-started")
# set(CPACK_NSIS_URL_INFO_ABOUT "https://www.openms.de")
# set(CPACK_NSIS_CONTACT "open-ms-general@lists.sourceforge.net")
# set(CPACK_NSIS_MENU_LINKS
#     "https://www.openms.de" "OpenMS Web Site")




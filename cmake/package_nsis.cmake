# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer, Chris Bielow $
# --------------------------------------------------------------------------


## Windows installer

# Check if nsis is actually run or if only the configuration is used
if("${CPACK_GENERATOR}" STREQUAL "NSIS")
  find_program(MAKENSIS_EXE makensis)

  if (NOT MAKENSIS_EXE)
    MESSAGE(FATAL_ERROR "Could not find 'makensis.exe'. Please make sure it's in $PATH!")
  endif()

  ## check for correct NSIS version
  execute_process(COMMAND ${MAKENSIS_EXE} /HDRINFO
                  OUTPUT_VARIABLE NSIS_INFO 
                  COMMAND_ERROR_IS_FATAL ANY)                
  STRING(FIND ${NSIS_INFO} "Size of each section is 16408 bytes" NSIS_IS_8K) ## the 1k version gives "2072 bytes"

  if (NSIS_IS_8K EQUAL -1)
    MESSAGE(FATAL_ERROR "NSIS (makensis.exe) needs to be the 'special build', which allows for 8k-length strings. This seems to be the 1k version. Please update NSIS. See https://github.com/OpenMS/NSIS")
  else()
    MESSAGE(STATUS "Found 8k version of NSIS. Great!")
  endif()
endif()
                

## With VS2019 the architecture HAS TO BE specified with the "–A" option or CMAKE_GENERATOR_PLATFORM var.
## Therefore the legacy way of adding a suffix to the Generator is not valid anymore.
## Read value of CMAKE_VS_PLATFORM_NAME instead
if (CMAKE_VS_PLATFORM_NAME MATCHES ".*Win32.*" OR CMAKE_GENERATOR MATCHES ".*Win32.*")
  set(PLATFORM "32")
  set(ARCH "x86")
else()
  set(PLATFORM "64")
  set(ARCH "x64")
endif()


#### Install System runtime libraries into /bin, so NSIS picks them up; this saves us from shipping a VC-Redist.exe with the installer
set(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE)
set (CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION ${INSTALL_LIB_DIR})
include(InstallRequiredSystemLibraries)


## Cfg_Settings.nsh only holds values that are fixed at configure time (redistributable
## paths, platform, source tree, version), so it can be written out directly here instead
## of being deferred to install time. It used to be generated in two stages: an @ONLY
## configure_file() here into a ".in.conf" file, then a *second* configure_file() inside
## an install(CODE ...) block (to substitute a PACKAGING_DIR only known once CPack starts
## staging files) that produced the final Cfg_Settings.nsh. That second stage depended on
## CPack's NSIS generator actually running our install(CODE) before it templated
## NSIS.template.in -- on Windows nightly builds it silently didn't (Cfg_Settings.nsh was
## never created at ${PROJECT_BINARY_DIR}, and "!include" failed with "could not find").
## The staging-directory-dependent OPENMS*DIR/THIRDPARTYDIR defines that used to live in
## Cfg_Settings.nsh.in are now defined directly in NSIS.template.in via CPack's own
## @CPACK_TEMPORARY_DIRECTORY@ substitution (the same mechanism already used for
## @CPACK_OPENMS_CFG_SETTINGS_FILE@ below), which CPack's NSIS generator fills in itself
## when it writes project.nsi -- no separate file generation step required.
configure_file(${PROJECT_SOURCE_DIR}/cmake/Windows/Cfg_Settings.nsh.in ${PROJECT_BINARY_DIR}/Cfg_Settings.nsh @ONLY)

## Pass the absolute path of the generated Cfg_Settings.nsh to the NSIS template via a CPACK_ variable
## (picked up automatically by include(CPack) below). We used to !include it via a hardcoded "..\..\..\"
## relative path from CPack's NSIS staging directory, which silently breaks whenever CPACK_PACKAGE_DIRECTORY
## is redirected elsewhere (e.g. to a short path on Windows to avoid MAX_PATH issues).
## Keep forward slashes here (like CPack's own CPACK_TOPLEVEL_DIRECTORY et al.): CPack re-serializes this
## variable into a quoted string in the generated CPackConfig.cmake, and a native Windows path's backslashes
## (e.g. "\a" in "D:\a\...", GitHub Actions' checkout root) are parsed as invalid CMake escape sequences there.
## NSIS itself accepts forward slashes in !include paths just fine.
set(CPACK_OPENMS_CFG_SETTINGS_FILE "${PROJECT_BINARY_DIR}/Cfg_Settings.nsh")

## Expose the install subdirectory layout to the NSIS template too (see NSIS.template.in,
## which combines these with CPack's own @CPACK_TEMPORARY_DIRECTORY@ to define the
## OPENMS*DIR/THIRDPARTYDIR paths). Only CPACK_-prefixed variables are persisted into
## CPackConfig.cmake and thus visible to CPack's template substitution at packaging time.
set(CPACK_OPENMS_INSTALL_BIN_DIR "${INSTALL_BIN_DIR}")
set(CPACK_OPENMS_INSTALL_LIB_DIR "${INSTALL_LIB_DIR}")
set(CPACK_OPENMS_INSTALL_SHARE_DIR "${INSTALL_SHARE_DIR}")
set(CPACK_OPENMS_INSTALL_PLUGIN_DIR "${INSTALL_PLUGIN_DIR}")
set(CPACK_OPENMS_INSTALL_DOC_DIR "${INSTALL_DOC_DIR}")

## Remove the next three lines if you use the NSIS autogeneration feature at some point!
## For now it makes sure everything is merged into the usual folders bin/share/include
set(CPACK_COMPONENT_ALL_IN_ONE 1)
set(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
set(CPACK_MONOLITHIC_INSTALL 1)
##

if((DEFINED ENV{CPACK_PACKAGE_FILE_NAME}) AND (NOT "$ENV{CPACK_PACKAGE_FILE_NAME}" STREQUAL ""))
  set(CPACK_PACKAGE_FILE_NAME "$ENV{CPACK_PACKAGE_FILE_NAME}")
else()
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}-Win${PLATFORM}")
endif()
set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS.ico")

## Create own target because you cannot "depend" on the internal target 'package'
add_custom_target(dist
  COMMAND ${CMAKE_COMMAND} -E echo "=== DIAGNOSTIC: Cfg_Settings.nsh state right before cpack invocation (inside dist target) ==="
  COMMAND cmd /c "if exist \"${PROJECT_BINARY_DIR}\\Cfg_Settings.nsh\" (echo FOUND ${PROJECT_BINARY_DIR}\\Cfg_Settings.nsh) else (echo MISSING ${PROJECT_BINARY_DIR}\\Cfg_Settings.nsh)"
  COMMAND cpack -G ${CPACK_GENERATOR} --verbose --debug
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




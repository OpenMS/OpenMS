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


## All NSIS.template.in settings that used to live in a separately !include'd
## Cfg_Settings.nsh are now defined directly in NSIS.template.in via CPack's own
## @CPACK_...@ template substitution (the same proven mechanism already used
## elsewhere in that file for e.g. @CPACK_TOPLEVEL_DIRECTORY@, @CPACK_PACKAGE_ICON@-
## derived paths, and the license page path -- all of which are absolute paths on
## the build tree's own drive and substitute reliably).
##
## Cfg_Settings.nsh was previously generated in two stages: an @ONLY configure_file()
## into a ".in.conf" file, then a second configure_file() deferred into an
## install(CODE ...) block (to substitute a PACKAGING_DIR only known once CPack
## starts staging files), and finally !include'd by NSIS.template.in via an
## absolute path. That worked right up until CPack actually ran makensis: CPack's
## own --debug trace confirmed the file existed, correctly populated, right up to
## (and including) the "if exist" check made microseconds before invoking cpack --
## yet NSIS's !include, running a fraction of a second later inside that same cpack
## invocation, reported it missing, with nothing in CPack's own instrumented file
## operations ever touching it. Whatever removes it, it is specific to a *separate*
## file surviving from CMake-configure-time through to NSIS-compile-time; values
## substituted directly into project.nsi by CPack's own template engine (as done
## here) never depend on any such file existing on disk at packaging time at all.
##
## Only CPACK_-prefixed variables are persisted into CPackConfig.cmake and thus
## visible to CPack's template substitution at packaging time.
set(CPACK_OPENMS_PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}")
set(CPACK_OPENMS_VERSION "${OPENMS_PACKAGE_VERSION_FULLSTRING}")
set(CPACK_OPENMS_PLATFORM "${PLATFORM}")
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




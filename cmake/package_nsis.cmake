# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------

## Windows installer

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
set(VC_REDIST_EXE "vcredist_${ARCH}.exe")


## Find redistributable to be installed by NSIS
if (NOT VC_REDIST_PATH)
	string(REGEX REPLACE ".*Visual Studio ([1-9][1-9]) .*" "\\1" OPENMS_MSVC_VERSION_STRING "${CMAKE_GENERATOR}")
	if("${OPENMS_MSVC_VERSION_STRING}" GREATER "14")
	  ## An old solution to find redists
	  #execute_process(COMMAND "$ENV{PROGRAMFILES}/Microsoft Visual Studio/Installer/vswhere" -latest -version "${OPENMS_MSVC_VERSION_STRING}" -property installationPath
	  #                OUTPUT_VARIABLE VC_ROOT_PATH
	  #				  ERROR_VARIABLE VSWHERE_ERROR
	  #				  RESULT_VARIABLE VSWHERE_RESULT)
	  #if ("${VSWHERE_RESULT}" NOTEQUAL "0")
	  #  message(FATAL_ERROR "Executing vswhere to find vsredist executable for win packaging failed. Either specify VC_REDIST_PATH or make sure vswhere works.")
	  #endif()
	  
	  ## We have to glob recurse in the parent folder because there is a version number in the end.
	  ## Unfortunately in my case the default version (latest) does not include the redist?!
	  ## TODO Not sure if this environment variable always exists. In the VS command line it should! Fallback vswhere or VCINSTALLDIR/Redist/MSVC?
	  get_filename_component(VC_ROOT_PATH "$ENV{VCToolsRedistDir}.." ABSOLUTE)
	  file(GLOB_RECURSE VC_REDIST_ABS_PATH "${VC_ROOT_PATH}/${VC_REDIST_EXE}")
	  ## TODO pick the latest of the found redists
	  get_filename_component(VC_REDIST_PATH "${VC_REDIST_ABS_PATH}" DIRECTORY)
	elseif(OPENMS_MSVC_VERSION_STRING GREATER "10")
	  set(VC_REDIST_PATH "$ENV{VCINSTALLDIR}redist/1033")  
	else()
	  message(FATAL_ERROR "Variable VC_REDIST_PATH missing."
	  "Before Visual Studio 2012 you have to provide the path"
	  "to the redistributable package of the VS you are using on your own.")
	endif()
endif()

if(EXISTS ${SEARCH_ENGINES_DIRECTORY})
  file(GLOB PWIZ_VCREDIST "${SEARCH_ENGINES_DIRECTORY}/*.exe")
  install(FILES ${PWIZ_VCREDIST}
          DESTINATION ${INSTALL_SHARE_DIR}/THIRDPARTY
		  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
		              GROUP_READ GROUP_EXECUTE
					  WORLD_READ WORLD_EXECUTE
		)
endif()

##TODO try following instead once CMake generates NSIS commands for us. Installs dll instead of redist though. Thirdparties?
# ########################################################### System runtime libraries
# set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP TRUE)
# include(InstallRequiredSystemLibraries)
# install(PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
#         DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/${PACKAGE_LIB_DIR}/
#         COMPONENT library)

## Careful: the configured file needs to lie exactly in the Build directory so that it is found by the NSIS_template
configure_file(${PROJECT_SOURCE_DIR}/cmake/Windows/Cfg_Settings.nsh.in ${PROJECT_BINARY_DIR}/Cfg_Settings.nsh.in.conf @ONLY)
install(CODE "
	set (PACKAGING_DIR \${CMAKE_INSTALL_PREFIX})
	configure_file(${PROJECT_BINARY_DIR}/Cfg_Settings.nsh.in.conf ${PROJECT_BINARY_DIR}/Cfg_Settings.nsh)
	")

set(CPACK_GENERATOR NSIS)
## Remove the next three lines if you use the NSIS autogeneration feature at some point!
## For now it makes sure everything is merged into the usual folders bin/share/include
set(CPACK_COMPONENT_ALL_IN_ONE 1)
set(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
set(CPACK_MONOLITHIC_INSTALL 1)
##

set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}-Win${PLATFORM}")
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




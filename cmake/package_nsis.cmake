# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2023.
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

## check if we are packaging at least Qt 5.15 (5.14 may also work but is untested), which is "-relocatable", i.e. can find ./bin/plugins/platforms/qwindows.dll without a qt.conf (which we do not ship anymore)
message(STATUS "Packaging: Checking Qt version ... found: ${Qt5Core_VERSION}")
if (Qt5Core_VERSION VERSION_LESS 5.15.0)
  message(FATAL_ERROR "Minimum supported Qt5 version is 5.15!")
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

if (NOT VC_REDIST_EXE)
	set(VC_REDIST_EXE "vcredist_${ARCH}.exe")
endif()

# Find redistributable to be installed by NSIS
if (NOT VC_REDIST_PATH)
	if (MSVC_TOOLSET_VERSION EQUAL 141)
		set(VS_VERSION 15)
	elseif(MSVC_TOOLSET_VERSION EQUAL 142)
		set(VS_VERSION 16)
	elseif(MSVC_TOOLSET_VERSION EQUAL 143)
		set(VS_VERSION 17)
	endif()

	if (DEFINED ENV{VCINSTALLDIR})
		## according to https://docs.microsoft.com/de-de/cpp/windows/redistributing-visual-cpp-files?view=msvc-160
		get_filename_component(VC_ROOT_PATH "$ENV{VCINSTALLDIR}Redist/MSVC/v${MSVC_TOOLSET_VERSION}" ABSOLUTE)
		file(GLOB_RECURSE VC_REDIST_ABS_PATH "${VC_ROOT_PATH}/vc_redist.${ARCH}.exe")
		if (NOT VC_REDIST_ABS_PATH)
			file(GLOB_RECURSE VC_REDIST_ABS_PATH "${VC_ROOT_PATH}/vcredist_${ARCH}.exe")
		endif()
		if (NOT VC_REDIST_ABS_PATH)
			file(GLOB_RECURSE VC_REDIST_ABS_PATH "${VC_ROOT_PATH}/vcredist.${ARCH}.exe")
		endif()

		message(STATUS "VS Redist locations (1st try): '${VC_REDIST_ABS_PATH}'")
	endif()

	if (NOT VC_REDIST_ABS_PATH AND DEFINED ENV{VCToolsRedistDir}) ## if still not found, try to use older methods
		## We have to glob recurse in the parent folder because there is a version number in the end.
		## Unfortunately in my case the default version (latest) does not include the redist?!
		## TODO Not sure if this environment variable always exists. In the VS command line it should! Fallback vswhere or VCINSTALLDIR/Redist/MSVC?
		get_filename_component(VC_ROOT_PATH "$ENV{VCToolsRedistDir}.." ABSOLUTE)
		file(GLOB_RECURSE VC_REDIST_ABS_PATH "${VC_ROOT_PATH}/vc_redist.${ARCH}.exe")
		if (NOT VC_REDIST_ABS_PATH)
			file(GLOB_RECURSE VC_REDIST_ABS_PATH "${VC_ROOT_PATH}/vcredist_${ARCH}.exe")
		endif()
		if (NOT VC_REDIST_ABS_PATH)
			file(GLOB_RECURSE VC_REDIST_ABS_PATH "${VC_ROOT_PATH}/vcredist.${ARCH}.exe")
		endif()

		message(STATUS "VS Redist locations (2nd try): '${VC_REDIST_ABS_PATH}'")
	endif()

	if (NOT VC_REDIST_ABS_PATH) ## if still not found try vswhere

		include(${OPENMS_HOST_DIRECTORY}/cmake/Windows/VSWhere.cmake)

		toolchain_ensure_vswhere()

		execute_process(COMMAND ${VSWHERE_PATH} -products * -latest -version "${VS_VERSION}" -property installationPath
		OUTPUT_VARIABLE VC_ROOT_PATH
		ERROR_VARIABLE VSWHERE_ERROR
		RESULT_VARIABLE VSWHERE_RESULT
		COMMAND_ECHO STDOUT
		)
		string(STRIP ${VC_ROOT_PATH} VC_ROOT_PATH)
		cmake_path(SET VC_ROOT_PATH NORMALIZE ${VC_ROOT_PATH})

		if (VSWHERE_RESULT EQUAL 0)
			message("Globbing for vc_redist.${ARCH}.exe in ${VC_ROOT_PATH}")
			file(GLOB_RECURSE VC_REDIST_ABS_PATH FOLLOW_SYMLINKS "${VC_ROOT_PATH}/vc_redist.${ARCH}.exe")
			if (NOT VC_REDIST_ABS_PATH)
				file(GLOB_RECURSE VC_REDIST_ABS_PATH "${VC_ROOT_PATH}/vcredist_${ARCH}.exe")
			endif()
			if (NOT VC_REDIST_ABS_PATH)
				file(GLOB_RECURSE VC_REDIST_ABS_PATH "${VC_ROOT_PATH}/vcredist.${ARCH}.exe")
			endif()
			message(STATUS "VS Redist locations (3rd try): ${VC_REDIST_ABS_PATH}")
		endif()
	endif()

	if (VC_REDIST_ABS_PATH)
		## arbitrarily pick first of the found redists
		list(GET VC_REDIST_ABS_PATH 0 VC_REDIST_ABS_PATH)
		get_filename_component(VC_REDIST_PATH "${VC_REDIST_ABS_PATH}" DIRECTORY)
		get_filename_component(VC_REDIST_EXE "${VC_REDIST_ABS_PATH}" NAME)
		message(STATUS "   ... picked first directory: ${VC_REDIST_PATH}")
	endif()
endif()

if (NOT VC_REDIST_PATH) ## if still not found: error
	message(FATAL_ERROR "Variable VC_REDIST_PATH missing. Are you on a proper VS 2019+ command line?")
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

#TODO try following instead once CMake generates NSIS commands for us. Installs dll instead of redist though.
# This means we would need to change the NSIS script to just copy them over.

# ########################################################### System runtime libraries
# Multiple ways to achieve that:
# - Currently: We package the vc_Redist installer exe from the VS directory
# - In package_general, we could adapt PRE and POST_EX/INCLUDES as well as maybe
#   DIRECTORIES to include them (from System32? or from VS directory?)
# - The code below also might achieve it.
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




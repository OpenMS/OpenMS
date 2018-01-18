# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
# $Authors: Stephan Aiche, Julianus Pfeuffer $
# --------------------------------------------------------------------------

# custom code to add a mac bundle application

# --------------------------------------------------------------------------
# adds a bundle with the given name
# @param _name Name of the bundle to add
# Note that this macro will also take care of installing the bundle in the
# case that a dmg package should be build
macro(add_mac_app_bundle _name)
	# the icon file
	set(ICON_FILE_PATH      "${PROJECT_SOURCE_DIR}/source/VISUAL/APPLICATIONS/GUITOOLS/${_name}-resources/${_name}.icns")
	set(INFO_PLIST_TEMPLATE "${PROJECT_SOURCE_DIR}/source/VISUAL/APPLICATIONS/GUITOOLS/${_name}-resources/${_name}.plist.in")
	get_filename_component(ICON_FILE_NAME "${ICON_FILE_PATH}" NAME)

	# we also need the icns in the app
	add_executable(
		${_name}
		MACOSX_BUNDLE
		${GUI_DIR}/${_name}.cpp
		${ICON_FILE_PATH})

	string(TIMESTAMP MY_YEAR "%Y")

	set_target_properties(${_name} PROPERTIES
		# we want our own info.plist template
		MACOSX_BUNDLE_INFO_PLIST "${INFO_PLIST_TEMPLATE}"
		MACOSX_BUNDLE_INFO_STRING "${PROJECT_NAME} Version ${CF_OPENMS_PACKAGE_VERSION}, Copyright ${MY_YEAR} The OpenMS Team."
		MACOSX_BUNDLE_ICON_FILE ${ICON_FILE_NAME}
		MACOSX_BUNDLE_GUI_IDENTIFIER "de.openms.${_name}"
		MACOSX_BUNDLE_LONG_VERSION_STRING "${PROJECT_NAME} Version ${CF_OPENMS_PACKAGE_VERSION}"
		MACOSX_BUNDLE_BUNDLE_NAME ${_name}
		MACOSX_BUNDLE_SHORT_VERSION_STRING ${CF_OPENMS_PACKAGE_VERSION}
		MACOSX_BUNDLE_BUNDLE_VERSION ${CF_OPENMS_PACKAGE_VERSION}
		MACOSX_BUNDLE_COPYRIGHT "Copyright ${MY_YEAR}, The OpenMS Team. All Rights Reserved."
	)

	set_source_files_properties(${ICON_FILE_PATH} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)

	## If you are packaging: Fix up the bundle. -> Copies all non-system dylibs into the bundles
	## Results in duplicate libraries in the different bundles.. but well.. that's how it is
	## If you are not packaging, libraries are linked via hardcoded paths specific to your machine.
	if("${PACKAGE_TYPE}" STREQUAL "dmg")
		install(CODE "
			set(BU_CHMOD_BUNDLE_ITEMS On)
			include(BundleUtilities)
			fixup_bundle(${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_name}.app \"\" \"\")
			"
			COMPONENT AApplications)

		install(TARGETS ${_name} BUNDLE
			DESTINATION .
			COMPONENT Applications)
		
		if(DEFINED CPACK_BUNDLE_APPLE_CERT_APP)
		   ## TODO try to find codesign to make sure the right exec is used
		   ## TODO allow choosing keychain
		   ## Note: Signing identity has to be unique, and present in any of the keychains in search list
		   ## which needs to be unlocked. Play around with keychain argument otherwise.
                   install(CODE "
execute_process(COMMAND codesign --deep --force --sign ${CPACK_BUNDLE_APPLE_CERT_APP} -i de.openms.${_name} \${CMAKE_INSTALL_PREFIX}/${_name}.app OUTPUT_VARIABLE sign_out ERROR_VARIABLE sign_out)
message('\${sign_out}')" COMPONENT BApplications)

                   install(CODE "
execute_process(COMMAND codesign -dv \${CMAKE_INSTALL_PREFIX}/${_name}.app OUTPUT_VARIABLE sign_check_out ERROR_VARIABLE sign_check_out)
message('\${sign_check_out}')" COMPONENT BApplications)
		 
		endif(DEFINED CPACK_BUNDLE_APPLE_CERT_APP)

	else()
	  ## Just install to the usual bin dir without fixing it up
	  install_tool(${_name})
	endif()
endmacro()

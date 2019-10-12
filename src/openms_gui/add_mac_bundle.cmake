# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
        set (APP_FOLDER "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_name}.app")

        set (PLUGIN_VAR_NAME QT_PLUGINS_APPS_${_name})
        ## Install qt5 plugins needed on mac and save them in a var for fixing their dependencies later
        install_qt5_plugin("Qt5::QCocoaIntegrationPlugin" ${PLUGIN_VAR_NAME} "${APP_FOLDER}/Contents/PlugIns" AAApplications)
        install_qt5_plugin("Qt5::QMacStylePlugin" ${PLUGIN_VAR_NAME} "${APP_FOLDER}/Contents/PlugIns" AAApplications)

        set (QT_PLUGINS_TO_FIX ${${PLUGIN_VAR_NAME}})

        ## Find Qt library folder
        get_target_property(QT_LIBRARY_DIR Qt5::Core LOCATION)
        get_filename_component(QT_LIBRARY_DIR ${QT_LIBRARY_DIR} PATH)
        get_filename_component(QT_LIBRARY_DIR "${QT_LIBRARY_DIR}/.." ABSOLUTE)

        ## Write a qt.conf file with a ref to the plugin dir in app bundles = PlugIns
        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
        "[Paths]\nPlugins = PlugIns\n")
        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
                DESTINATION "${APP_FOLDER}/Contents/Resources"
                COMPONENT AAApplications)

        ## Fix up the dependencies in the bundle and make them rel. to their location in the bundle
        ## Give additional plugins to fix and extra dirs where dependencies should be searched
		install(CODE "
			set(BU_CHMOD_BUNDLE_ITEMS On)
			include(BundleUtilities)
			fixup_bundle(${APP_FOLDER} \"${QT_PLUGINS_TO_FIX}\" \"${QT_LIBRARY_DIR}\")
			"
			COMPONENT AApplications)

		## Copy bundle to the target install destination
		install(TARGETS ${_name} BUNDLE
			DESTINATION .
			COMPONENT Applications)
		
		## Notarization is only possible with Xcode/Appleclang 10, otherwise we just skip
		if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
			## We also need an identity to sign with
			if(DEFINED CPACK_BUNDLE_APPLE_CERT_APP)
			## TODO try to find codesign to make sure the right exec is used (currently needs to be in path)
			## TODO allow choosing keychain
			## Note: Signing identity has to be unique, and present in any of the keychains in search list
			## which needs to be unlocked. Play around with keychain argument otherwise.
                   install(CODE "
execute_process(COMMAND codesign --deep --force --sign ${CPACK_BUNDLE_APPLE_CERT_APP} -i de.openms.${_name} \${CMAKE_INSTALL_PREFIX}/${_name}.app OUTPUT_VARIABLE sign_out ERROR_VARIABLE sign_out)
message('\${sign_out}')" COMPONENT BApplications)

                   install(CODE "
execute_process(COMMAND codesign -dv \${CMAKE_INSTALL_PREFIX}/${_name}.app OUTPUT_VARIABLE sign_check_out ERROR_VARIABLE sign_check_out)
message('\${sign_check_out}')" COMPONENT BApplications)

                   install(CODE "
execute_process(COMMAND ditto -c -k --rsrc --keepParent \${CMAKE_INSTALL_PREFIX}/${_name}.app \${CMAKE_INSTALL_PREFIX}/${_name}.app.zip OUTPUT_VARIABLE ditto_out ERROR_VARIABLE ditto_out)
message('\${ditto_out}')" COMPONENT BApplications)

                   install(CODE "
execute_process(COMMAND ${OPENMS_HOST_DIRECTORY}/cmake/MacOSX/notarize_app.sh \${CMAKE_INSTALL_PREFIX}/${_name}.app.zip de.openms.${_name} ${SIGNING_EMAIL} CODESIGNPW ${OPENMS_HOST_BINARY_DIRECTORY} OUTPUT_VARIABLE notarize_out ERROR_VARIABLE notarize_out)
message('\${notarize_out}')" COMPONENT BApplications)

                   install(CODE "
execute_process(COMMAND spctl -a -v \${CMAKE_INSTALL_PREFIX}/${_name}.app OUTPUT_VARIABLE verify_out ERROR_VARIABLE verify_out)
message('\${verify_out}')" COMPONENT BApplications)

                    install(CODE "
execute_process(COMMAND rm \${CMAKE_INSTALL_PREFIX}/${_name}.app.zip OUTPUT_VARIABLE remove_out ERROR_VARIABLE remove_out)
message('\${remove_out}')" COMPONENT BApplications)
                   
			endif(DEFINED CPACK_BUNDLE_APPLE_CERT_APP)
		endif()
	else()
	  ## Just install to the usual bin dir without fixing it up
	  install_tool(${_name})
	endif()
endmacro()

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
# $Authors: Stephan Aiche, Julianus Pfeuffer $
# --------------------------------------------------------------------------

## Very useful for debugging purposes: Disables the dependency on the "make install" target.
## In our case e.g. "make install" always builds the documentation etc. 
#set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY On)
## Or reduce the set of components to package
#set(CPACK_COMPONENTS_ALL Applications Dependencies zzz-fixing-dependencies)


## We want to package the whole top-level dir so a user can drag'n'drop it via the image.
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 0) ## dmg seems to be component-aware and makes an ALL-IN-ONE package
set(CPACK_COMPONENT_INCLUDE_TOPLEVEL_DIRECTORY 1) ## Therefore _only_ use the second.. weird stuff.

set(MACOS_TARGET_ARCHS ${CMAKE_OSX_ARCHITECTURES})
if (NOT MACOS_TARGET_ARCHS)
  # Warning: if cmake is a subprocess of a process that is run under Rosetta, it will yield
  #  x86_64 (but probably also build for it. Therefore it should be fine.)
  set(MACOS_TARGET_ARCHS ${CMAKE_HOST_SYSTEM_PROCESSOR})
endif()
if (MACOS_TARGET_ARCHS STREQUAL "x86_64")
  set(ARCH_SUFFIX "Intel")
elseif (MACOS_TARGET_ARCHS STREQUAL "arm64")
  set(ARCH_SUFFIX "Silicon")
elseif ("x86_64" IN_LIST MACOS_TARGET_ARCHS AND "arm64" IN_LIST MACOS_TARGET_ARCHS)
  set(ARCH_SUFFIX "Universal")
else ()
  set(ARCH_SUFFIX "Unknown")
  message(WARNING "Couldn't determine MACOS_TARGET_ARCHS.")
endif()

if((DEFINED ENV{CPACK_PACKAGE_FILE_NAME}) AND (NOT "$ENV{CPACK_PACKAGE_FILE_NAME}" STREQUAL ""))
  set(CPACK_PACKAGE_FILE_NAME "$ENV{CPACK_PACKAGE_FILE_NAME}")
else()
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}-macOS-${ARCH_SUFFIX}")
endif()

## Note: That the mac app bundles (TOPPView) take care of themselves
##       when installed as dmg (see src/openms_gui/add_mac_bundle.cmake)

## Fix OpenMS dependencies for all executables in the install directory under bin.
## That affects everything but the bundles (which have their own structure and are fixed up in add_mac_bundle.cmake)
########################################################### Fix Dependencies
install(CODE "execute_process(COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_dependencies.rb -b \${CMAKE_INSTALL_PREFIX}/${INSTALL_BIN_DIR}/ -l \${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/)"
  COMPONENT zzz-fixing-dependencies
) # not applicable as long as we dont use Qt plugins in the TOPP tools: -p \${CMAKE_INSTALL_PREFIX}/${INSTALL_PLUGIN_DIR}) 

if(DEFINED CPACK_BUNDLE_APPLE_CERT_APP)
  install(CODE "execute_process(COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/sign_bins_and_libs.rb -d \${CMAKE_INSTALL_PREFIX}/ -s ${CPACK_BUNDLE_APPLE_CERT_APP})"
    COMPONENT zzz-sign-bins-and-libs
  )
endif()

## Additionally install TOPPShell into root of install folder

########################################################### TOPPShell
install(FILES       ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP-shell.command
        DESTINATION .
        RENAME      "TOPP Shell"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_READ GROUP_EXECUTE
                    WORLD_READ WORLD_EXECUTE
        COMPONENT   TOPPShell)

install(FILES       ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP_bash_profile
        DESTINATION .
        RENAME      .TOPP_bash_profile
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
        COMPONENT   TOPPShell)

install(FILES       ${PROJECT_SOURCE_DIR}/cmake/MacOSX/README.md
        DESTINATION .
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
        COMPONENT   TOPPShell)

## ----------
## No need to for this at the moment: no Qt plugins for the general CLI tools required anymore. Left here for reference.
## ----------
## Install the qt.conf file so we can find the libraries
## add qt.conf to the bin directory for DMGs
#file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
#"[Paths]\nPlugins = ../${INSTALL_PLUGIN_DIR}/\n")
#install(FILES       ${CMAKE_CURRENT_BINARY_DIR}/qt.conf
#        DESTINATION ${INSTALL_BIN_DIR}
#        PERMISSIONS OWNER_WRITE OWNER_READ
#                    GROUP_READ
#                    WORLD_READ
#        COMPONENT   Dependencies)


## Create own target because you cannot "depend" on the internal target 'package'
add_custom_target(dist
  COMMAND cpack -G ${CPACK_GENERATOR}
  COMMENT "Building ${CPACK_GENERATOR} package"
)

########################################################### Create dmg with background image
if (DEFINED CMAKE_VERSION AND NOT "${CMAKE_VERSION}" VERSION_LESS "3.5")
  set(OPENMS_DMG_FOLDER_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}") ## The name of the OpenMS folder on the DMG
  configure_file(${PROJECT_SOURCE_DIR}/cmake/MacOSX/setup_applescript.scpt.in ${PROJECT_BINARY_DIR}/macOS_bundle_setup/setup_applescript.scpt)
  set(CPACK_DMG_DS_STORE_SETUP_SCRIPT ${PROJECT_BINARY_DIR}/macOS_bundle_setup/setup_applescript.scpt)
  #Next line could overcome a script but since we do not have a fixed name of the OpenMS-$VERSION folder, it probably won't work
  #set(CPACK_DMG_DS_STORE ${PROJECT_SOURCE_DIR}/cmake/MacOSX/DS_store_new)
  set(CPACK_DMG_BACKGROUND_IMAGE ${PROJECT_SOURCE_DIR}/cmake/MacOSX/background.png)
  set(CPACK_DMG_FORMAT UDBZ) ## Try bzip2 to get slightly smaller images
  
  ## Sign the image. CPACK_BUNDLE_APPLE_CERT_APP needs to be unique and found in one of the
  ## keychains in the search list (which needs to be unlocked).
  if (DEFINED CPACK_BUNDLE_APPLE_CERT_APP)
    add_custom_target(signed_dist
                      COMMAND codesign --deep --force --sign ${CPACK_BUNDLE_APPLE_CERT_APP} ${CPACK_PACKAGE_FILE_NAME}.dmg
                      COMMAND ${OPENMS_HOST_DIRECTORY}/cmake/MacOSX/notarize_app.sh ${CPACK_PACKAGE_FILE_NAME}.dmg de.openms ${SIGNING_EMAIL} CODESIGNPW ${OPENMS_HOST_BINARY_DIRECTORY}
                      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                      COMMENT "Signing and notarizing ${CPACK_PACKAGE_FILE_NAME}.dmg as ${CPACK_BUNDLE_APPLE_CERT_APP}"
                      DEPENDS dist)
  endif()
  
else()
  ## The old scripts need the background image in the target folder.
  ########################################################### Background Image
  install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/background.png
        DESTINATION share/OpenMS/
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
        COMPONENT share)
  
  ## Next command assumes the dmg was already generated and lies in the build directory.
  add_custom_target(finalized_dist
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fixdmg.sh
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Finalizing dmg image"
    DEPENDS dist)
endif()

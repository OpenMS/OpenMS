# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

## Very useful for debugging purposes: Disables the dependency on the "make install" target.
## In our case e.g. "make install" always builds the documentation. 
#set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY On)

set(CPACK_GENERATOR "DragNDrop")

## drag'n'drop installation configuration
## Note: We have certain dependencies between the individual components!!!
##       To ensure that the components are executed in the correct order
##       we use the fact that cmake executes them in alphabetical order
##        1. A-Z
##        2. a-z
##       So before adding an additional target make sure that you do not
##       intefer with other namings/components.
##
## Note: That the mac app bundles (TOPPView) take care of them selfes
##       when installed as dmg (see src/openms_gui/add_mac_bundle.cmake)

## If you want to rename the binary dir on the dmg.
set(DMG_BINARY_DIR_NAME "TOPP")

## Fix OpenMS dependencies for all executables in the install directory
########################################################### Fix Dependencies
install(CODE "execute_process(COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_dependencies.rb -b \${CMAKE_INSTALL_PREFIX}/${CPACK_PACKAGE_INSTALL_DIRECTORY}/${DMG_BINARY_DIR_NAME}/ -l \${CMAKE_INSTALL_PREFIX}/${CPACK_PACKAGE_INSTALL_DIRECTORY}/lib/ -v)"
  COMPONENT zzz-fixing-dependencies
)

########################################################### Libraries
# Libraries hack, avoid cmake interferring with our own lib fixes
install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/
  DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/lib/
  COMPONENT library
)

########################################################### TOPP Binaries
# Binary hack, avoid cmake interferring with our own lib fixes
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin/
	DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/${DMG_BINARY_DIR_NAME}
  COMPONENT applications
  FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_READ GROUP_EXECUTE
                        WORLD_READ WORLD_EXECUTE
  DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_READ GROUP_EXECUTE
                        WORLD_READ WORLD_EXECUTE
  PATTERN "*.app" EXCLUDE
  REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
  REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
)

########################################################### Share
install(DIRECTORY share/
	DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/share
	COMPONENT share
  FILE_PERMISSIONS      OWNER_WRITE OWNER_READ
                        GROUP_READ
                        WORLD_READ
  DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_EXECUTE GROUP_READ
                        WORLD_EXECUTE WORLD_READ
  REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
  REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
)

########################################################### Documentation Preparation
#------------------------------------------------------------------------------
# Since doc_tutorial is built with ALL when ENABLE_TUTORIALS is ON, we can assume these
# PDFs are present. The ALL target is executed automatically with "make package"

if (ENABLE_TUTORIALS)
  install(FILES       ${PROJECT_BINARY_DIR}/doc/OpenMS_tutorial.pdf
          DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/
          COMPONENT doc
          PERMISSIONS OWNER_WRITE OWNER_READ
                      GROUP_READ
                      WORLD_READ

  )

  install(FILES     ${PROJECT_BINARY_DIR}/doc/TOPP_tutorial.pdf
          DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/
          COMPONENT doc
          PERMISSIONS OWNER_WRITE OWNER_READ
                      GROUP_READ
                      WORLD_READ
  )
else()
  message("Warning: Configuring for packaging without tutorials. If you want to build a full package make sure that configuration with -DENABLE_TUTORIALS=On succeeds.")
endif ()


########################################################### Documentation
install(FILES       ${PROJECT_BINARY_DIR}/doc/index.html
        DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/
        RENAME OpenMSAndTOPPDocumentation.html
        COMPONENT doc
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
)

install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html
        DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/
        COMPONENT doc
        FILE_PERMISSIONS      OWNER_WRITE OWNER_READ
                              GROUP_READ
                              WORLD_READ
        DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                              GROUP_EXECUTE GROUP_READ
                              WORLD_EXECUTE WORLD_READ
        REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
        REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
)

########################################################### SEARCHENGINES
if(EXISTS ${SEARCH_ENGINES_DIRECTORY})
  if(EXISTS ${SEARCH_ENGINES_DIRECTORY}/Fido)
    install(DIRECTORY             ${SEARCH_ENGINES_DIRECTORY}/Fido
            DESTINATION           OpenMS-${CPACK_PACKAGE_VERSION}/TOPP/SEARCHENGINES
            COMPONENT             SearchEngine-Fido
            FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
            REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
            )
  endif()

  if(EXISTS ${SEARCH_ENGINES_DIRECTORY}/MSGFPlus)
    install(DIRECTORY             ${SEARCH_ENGINES_DIRECTORY}/MSGFPlus
            DESTINATION           OpenMS-${CPACK_PACKAGE_VERSION}/TOPP/SEARCHENGINES
            COMPONENT             SearchEngine-MSGFPlus
            FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
            REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
            )
  endif()

  if(EXISTS ${SEARCH_ENGINES_DIRECTORY}/OMSSA)
    install(DIRECTORY             ${SEARCH_ENGINES_DIRECTORY}/OMSSA
            DESTINATION           OpenMS-${CPACK_PACKAGE_VERSION}/TOPP/SEARCHENGINES
            COMPONENT             SearchEngine-OMSSA
            FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
            REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
            )
  endif()

  if(EXISTS ${SEARCH_ENGINES_DIRECTORY}/XTandem)
    install(DIRECTORY             ${SEARCH_ENGINES_DIRECTORY}/XTandem
            DESTINATION           OpenMS-${CPACK_PACKAGE_VERSION}/TOPP/SEARCHENGINES
            COMPONENT             SearchEngine-XTandem
            FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
            REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
            )
  endif()
  if(EXISTS ${SEARCH_ENGINES_DIRECTORY}/LuciPHOr2)
    install(DIRECTORY             ${SEARCH_ENGINES_DIRECTORY}/LuciPHOr2
            DESTINATION           OpenMS-${CPACK_PACKAGE_VERSION}/TOPP/SEARCHENGINES
            COMPONENT             SearchEngine-LuciPHOr2
            FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
            REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
            )
  endif()
  ## MyriMatch does not exist for MacOSX
endif()

########################################################### TOPPShell
install(FILES       ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP-shell.command
        DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/
        RENAME      "TOPP Shell"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_READ GROUP_EXECUTE
                    WORLD_READ WORLD_EXECUTE
        COMPONENT   TOPPShell)

install(FILES       ${PROJECT_SOURCE_DIR}/cmake/MacOSX/TOPP_bash_profile
        DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/
        RENAME      .TOPP_bash_profile
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
        COMPONENT   TOPPShell)

install(FILES       ${PROJECT_SOURCE_DIR}/cmake/MacOSX/README
        DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
        COMPONENT   TOPPShell)

########################################################### Create dmg with background image
if (DEFINED CMAKE_VERSION AND NOT "${CMAKE_VERSION}" VERSION_LESS "3.5")
  set(OPENMS_LOGO ${PROJECT_SOURCE_DIR}/cmake/MacOSX/openms_logo_large_transparent.png) ## For configuration of the script
  configure_file(${PROJECT_SOURCE_DIR}/cmake/MacOSX/setup_applescript.scpt.in ${PROJECT_BINARY_DIR}/macOS_bundle_setup/setup_applescript.scpt)
  set(CPACK_DMG_DS_STORE_SETUP_SCRIPT ${PROJECT_BINARY_DIR}/macOS_bundle_setup/setup_applescript.scpt)
  #Next line could overcome a script but since we do not have a fixed name of the OpenMS-$VERSION folder, it probably won't work
  #set(CPACK_DMG_DS_STORE ${PROJECT_SOURCE_DIR}/cmake/MacOSX/DS_store_new)
  set(CPACK_DMG_BACKGROUND_IMAGE ${PROJECT_SOURCE_DIR}/cmake/MacOSX/background.png)
else()
  ## The old scripts need the background image in the target folder.
  ########################################################### Background Image
  install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/background.png
        DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/share/OpenMS/
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
        COMPONENT share)
  ## Use custom scripts and targets (just 'make package' creates a very simple dmg).
  add_custom_target(final_package
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fixdmg.sh
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Finalizing dmg image"
    DEPENDS dmg)

  add_custom_target(dmg
    COMMAND cpack -G DragNDrop
    COMMENT "Building intermediate dmg package"
  )
endif()

include(CPack)


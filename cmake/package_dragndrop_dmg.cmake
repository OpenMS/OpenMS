# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

set(CPACK_GENERATOR "DragNDrop")

## drag'n'drop installaltion configuration
## Note: We have certain dependencies between the individual components!!!
##       To ensure that the components are executed in the correct order 
##       we use the fact that cmake executes them in alphabetical order 
##        1. A-Z 
##        2. a-z
##       So before adding an additional target make sure that you do not
##       intefer with other namings/components.

## Fix OpenMS dependencies for all executables
########################################################### Fix Dependencies
install(CODE "execute_process(COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_dependencies.rb -b ${PROJECT_BINARY_DIR}/bin/ -l ${PROJECT_BINARY_DIR}/lib/ -v)"
  COMPONENT Fixing-dependencies
)

# Install the apps
########################################################### TOPPAS, TOPPView, INIFileEditor
set(APP_BUNDLES "TOPPAS" "TOPPView" "INIFileEditor")

foreach(APP_BUNDLE IN LISTS APP_BUNDLES)
  install(CODE "
  	set(BU_CHMOD_BUNDLE_ITEMS On)
  	include(BundleUtilities)
  	fixup_bundle(${PROJECT_BINARY_DIR}/${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${APP_BUNDLE}.app \"\" \"\")"
  	COMPONENT AApplications)

  install(TARGETS ${APP_BUNDLE} BUNDLE
  				DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/
  				COMPONENT Applications)
endforeach()

########################################################### Libraries
# Libraries hack, avoid cmake interferring with our own lib fixes
install(DIRECTORY ${PROJECT_BINARY_DIR}/lib/
  DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/lib/
  COMPONENT library
)

########################################################### TOPP Binaries
# Binary hack, avoid cmake interferring with our own lib fixes
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin/
	DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/TOPP
  COMPONENT applications
  FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_READ GROUP_EXECUTE
                        WORLD_READ WORLD_EXECUTE
  DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_READ GROUP_EXECUTE
                        WORLD_READ WORLD_EXECUTE
  PATTERN "*.app" EXCLUDE
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
	PATTERN ".svn" EXCLUDE
)

########################################################### Documentation Preparation
install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build \"${PROJECT_BINARY_DIR}\" --target doc)"
        COMPONENT AAA-Documentation-Preparation)
install(CODE "execute_process(COMMAND \"${CMAKE_COMMAND}\" --build \"${PROJECT_BINARY_DIR}\" --target doc_tutorials)"
        COMPONENT AAA-Documentation-Preparation)
        
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
        PATTERN ".svn" EXCLUDE
)

install(FILES 		  ${PROJECT_BINARY_DIR}/doc/OpenMS_tutorial.pdf
        DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/
        COMPONENT doc
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
)

install(FILES 		${PROJECT_BINARY_DIR}/doc/TOPP_tutorial.pdf
        DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/Documentation/
        COMPONENT doc
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
)

########################################################### SEARCHENGINES
if(EXISTS ${SEARCH_ENGINES_DIRECTORY})
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
            )
  endif()
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

########################################################### Background Image
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MacOSX/background.png
        DESTINATION OpenMS-${CPACK_PACKAGE_VERSION}/share/OpenMS/
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
        COMPONENT share)

include(CPack)

########################################################### Create dmg with background image
add_custom_target(final_package
  COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fixdmg.sh
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  COMMENT "Finalizing dmg image"
  DEPENDS dmg)

add_custom_target(dmg
  COMMAND cpack -G DragNDrop
  COMMENT "Building intermediate dmg package"
)

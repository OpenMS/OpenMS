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

## Very useful for debugging purposes: Disables the dependency on the "make install" target.
## In our case e.g. "make install" always builds the documentation. 
#set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY On)

set(CPACK_GENERATOR "DragNDrop")

## We want to package the whole top-level dir so a user can drag'n'drop it via the image.
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 0) ## dmg seems to be component-aware and makes an ALL-IN-ONE package
set(CPACK_COMPONENT_INCLUDE_TOPLEVEL_DIRECTORY 1) ## Therefore _only_ use the second.. weird stuff.

## Note: That the mac app bundles (TOPPView) take care of themselves
##       when installed as dmg (see src/openms_gui/add_mac_bundle.cmake)


## Fix OpenMS dependencies for all executables in the install directory
########################################################### Fix Dependencies
install(CODE "execute_process(COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_dependencies.rb -b \${CMAKE_INSTALL_PREFIX}/${INSTALL_BIN_DIR}/ -l \${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/ -v)"
  COMPONENT zzz-fixing-dependencies
)

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

install(FILES       ${PROJECT_SOURCE_DIR}/cmake/MacOSX/README
        DESTINATION .
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
        DESTINATION share/OpenMS/
        PERMISSIONS OWNER_WRITE OWNER_READ
                    GROUP_READ
                    WORLD_READ
        COMPONENT share)
  ## Use custom scripts and targets (just 'make package' creates a very simple dmg).
  add_custom_target(dmg
    COMMAND cpack -G DragNDrop
    COMMENT "Building intermediate dmg package"
  )
  
  ## Next command assumes the dmg was already generated and lies in the build directory.
  add_custom_target(final_package
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fixdmg.sh
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Finalizing dmg image"
    DEPENDS dmg)
endif()

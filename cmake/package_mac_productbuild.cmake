# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright OpenMS Inc. -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-present.
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

## Very useful for debugging purposes: Disables the dependency on the "make install" target.
## In our case e.g. "make install" always builds the documentation etc. 
#set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY On)

# For info about this CPack generator and its capabilities, see the CMake documentation
# For info about productbuild and the flat package format see https://matthew-brett.github.io/docosx/flat_packages.html

set(CPACK_PACKAGING_INSTALL_PREFIX "/Applications/${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}")
set(CPACK_PRODUCTBUILD_IDENTIFIER "de.openms")
set(CPACK_PRODUCTBUILD_RESOURCES_DIR ${PROJECT_SOURCE_DIR}/cmake/MacOSX)
set(CPACK_PRODUCTBUILD_BACKGROUND ${OPENMS_LOGOSMALL_NAME})
set(CPACK_PRODUCTBUILD_BACKGROUND_ALIGNMENT "bottomleft")
set(CPACK_PRODUCTBUILD_BACKGROUND_SCALING "none")

# Allow installing to every Domain if supported by current CMake version (https://gitlab.kitware.com/cmake/cmake/-/merge_requests/6825)
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.23.0")
  set(CPACK_PRODUCTBUILD_DOMAINS TRUE) # system-wide
  set(CPACK_PRODUCTBUILD_DOMAINS_USER TRUE) # user folder
endif()

# TODO we might need to set a user-defined template for the installer anyway due to missing architecture support
# in CMake (https://gitlab.kitware.com/cmake/cmake/-/issues/21734)
# The template would go in cmake/Modules which is already in our Module path.
# Official template is here: https://gitlab.kitware.com/cmake/cmake/-/blob/v3.27.4/Modules/Internal/CPack/CPack.distribution.dist.in?ref_type=tags

if(NOT DEFINED CPACK_PRODUCTBUILD_IDENTITY_NAME)
  message(WARNING "CPACK_PRODUCTBUILD_IDENTITY_NAME not set. PKG will not be signed. Make sure to specify an identity with a Developer ID: Installer certificate (not Application certificate).")
endif()

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

## Not needed unless we need Qt plugins for TOPP again
## Install the qt.conf file so we can find the libraries
## add qt.conf to the bin directory for DMGs/pkgs
#file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/qt.conf"
#"[Paths]\nPlugins = ../${INSTALL_PLUGIN_DIR}\n")
#install(FILES       ${CMAKE_CURRENT_BINARY_DIR}/qt.conf
#        DESTINATION ./${INSTALL_BIN_DIR}
#        PERMISSIONS OWNER_WRITE OWNER_READ
#                    GROUP_READ
#                    WORLD_READ
#        COMPONENT   Applications)

## Fix OpenMS dependencies for all executables in the install directory under bin.
## That affects everything but the bundles (whose Framework folders are symlinked to lib anyway).
########################################################### Fix Dependencies
# Use --no-copy and assume CMake correctly collected all necessary libs already (therefore only fix links to dylibs)
install(CODE "execute_process(COMMAND ${OPENMS_HOST_DIRECTORY}/cmake/MacOSX/fix_dependencies.rb -b \${CMAKE_INSTALL_PREFIX}/../../../Applications${CPACK_PACKAGING_INSTALL_PREFIX}/${INSTALL_BIN_DIR}/ -l \${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/ -p \${CMAKE_INSTALL_PREFIX}/${INSTALL_PLUGIN_DIR} -e @rpath/ -n -c)"
        COMPONENT Dependencies
        )
# We have to do signing after fixing, otherwise signature will be invalidated.
# We have to choose COMPONENT Dependencies, otherwise it would be performed at the end of installing Applications, which
# comes first (alphabetically per component, then order of install calls)
# If the install CODE is not in the same component though, we need to navigate from the component specific install
# prefix to the other prefix. This is unfortunately very unrobust.
# TODO find better order or rewrite fix_dependencies script to be called separately
install(CODE "
        execute_process(COMMAND find \${CMAKE_INSTALL_PREFIX}/../../../Applications${CPACK_PACKAGING_INSTALL_PREFIX}/${INSTALL_BIN_DIR}/ -type f -execdir codesign --force --options runtime -i de.openms.TOPP.{} --sign \"${CPACK_BUNDLE_APPLE_CERT_APP}\" {} \\; OUTPUT_VARIABLE topp_sign_out ERROR_VARIABLE topp_sign_out)
        execute_process(COMMAND find \${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/ -type f -execdir codesign --force --options runtime -i de.openms.TOPP.libs.{} --sign \"${CPACK_BUNDLE_APPLE_CERT_APP}\" {} \\; OUTPUT_VARIABLE topp_sign_out ERROR_VARIABLE topp_sign_out)
        message('\${topp_sign_out}')"
        COMPONENT Dependencies
        )
install(CODE "execute_process(COMMAND ${OPENMS_HOST_DIRECTORY}/cmake/MacOSX/fix_dependencies.rb -l \${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/ -e @rpath/ -n -c)"
        COMPONENT library
        )
install(CODE "
        execute_process(COMMAND find \${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}/ -type f -execdir codesign --force --options runtime -i de.openms.TOPP.libs.{} --sign \"${CPACK_BUNDLE_APPLE_CERT_APP}\" {} \\; OUTPUT_VARIABLE lib_sign_out ERROR_VARIABLE lib_sign_out)
        message('\${lib_sign_out}')"
        COMPONENT library
        )

## When Applications are installed (which is the FIRST in alphabetical order AND the main component),
## a postinstall script runs to set file icon
install(FILES       ${PROJECT_SOURCE_DIR}/cmake/MacOSX/openms_logo_large_transparent.png
        DESTINATION .
        COMPONENT   Applications)
configure_file(${OPENMS_HOST_DIRECTORY}/cmake/MacOSX/setIcon.sh.in ${OPENMS_HOST_BINARY_DIRECTORY}/cmake/MacOSX/setIcon.sh)
set(CPACK_POSTFLIGHT_APPLICATIONS_SCRIPT ${OPENMS_HOST_BINARY_DIRECTORY}/cmake/MacOSX/setIcon.sh)


## Create own target because you cannot "depend" on the internal target 'package'
add_custom_target(dist
  COMMAND cpack -G ${CPACK_GENERATOR}
  COMMENT "Building ${CPACK_GENERATOR} package"
)

add_custom_target(finalized_dist
                  COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fixdmg.sh
                  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                  COMMENT "Finalizing dmg image"
                  DEPENDS dist)



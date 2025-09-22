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
#    documentation and/o>r other materials provided with the distribution.
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


# --------------------------------------------------------------------------
# general definitions used for building OpenMS packages
set(CPACK_PACKAGE_NAME "OpenMS")
set(CPACK_PACKAGE_VENDOR "OpenMS.de")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMS - A framework for mass spectrometry")
set(CPACK_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${OPENMS_PACKAGE_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${OPENMS_PACKAGE_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${OPENMS_PACKAGE_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageDescriptionFile.txt)
set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/License.txt)
set(CPACK_RESOURCE_FILE_WELCOME ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceWelcomeFile.txt)
set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceReadme.txt)
set(CPACK_STRIP_FILES TRUE) # to save some space in the installers

set(OPENMS_LOGO_NAME openms_logo_large_transparent.png) ## The filename of the logo to be used for the OpenMS folder e.g. on the DMG
set(OPENMS_LOGO ${PROJECT_SOURCE_DIR}/cmake/MacOSX/${OPENMS_LOGO_NAME}) ## The logo to be used for the OpenMS folder e.g. on the DMG

set(OPENMS_LOGOSMALL_NAME openms_logo_corner_transparent.png) ## The filename of the logo to be used for the OpenMS folder e.g. on the PKG
set(OPENMS_LOGOSMALL ${PROJECT_SOURCE_DIR}/cmake/MacOSX/${OPENMS_LOGOSMALL_NAME}) ## The logo to be used for the OpenMS folder e.g. on the PKG

########################################################### Fixing dynamic dependencies
## Qt Plugins needed for the CL tools should have been installed before (such as QSqliteDriverPlugin)
## This currently works because our libs and TOPP tools include all dependencies. For macOS,
##  the app bundles need to have a different RUNTIME_DEPENDENCY_SET (TOPPView_DEPS, ...) due
##  to CMake assuming you want standalone bundles. But we want to share libs between them.

# This would be to look in the Contrib and other cmake_prefix_paths for dependencies.
#list(TRANSFORM CMAKE_PREFIX_PATH APPEND "/bin" OUTPUT_VARIABLE DEP_BIN_DIRS)
#list(TRANSFORM CMAKE_PREFIX_PATH APPEND "/lib" OUTPUT_VARIABLE DEP_LIB_DIRS)
# But since we copy them in the build stage to our runtime directory (bin), we can add this one.

# On Windows we need to tell CMake where to look for.
# We also do not need API sets. So exclude them.
# Define shared exclusion patterns for all RUNTIME_DEPENDENCY_SET usage
if(WIN32)
  # Skip system runtime libraries detection
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP TRUE)
  
  # Enhanced exclusion patterns for Windows security and system DLLs
  set(OPENMS_RUNTIME_EXCLUDE "api-ms" "ext-ms" "hvsi" "pdmutilities" "wpaxholder" "dxgi" "uxtheme" "d3d11" "winmm" "wldp"
                              "AzureAttestManager" "AzureAttestNormal" "WTDSENSOR" "wtdccm" "bcrypt" "ntdll" "kernel32"
                              "advapi32" "user32" "ws2_32" "secur32" "crypt32" "wintrust" "msvcrt")
  set(OPENMS_RUNTIME_POST_EXCLUDE ".*[Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\\\]system32.*"
                                  ".*[Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\\\][Ss]ys[Ww][Oo][Ww]64.*"
                                  ".*[Aa]zure[Aa]ttest.*" ".*[Ww][Tt][Dd].*" ".*bcrypt.*" ".*ntdll.*"
                                  ".*kernel32.*" ".*advapi32.*" ".*user32.*" ".*ws2_32.*" ".*secur32.*"
                                  ".*crypt32.*" ".*wintrust.*" ".*msvcrt.*")
  # Backward compatibility
  set(EXCLUDE ${OPENMS_RUNTIME_EXCLUDE})
  set(POST_EXCLUDE ${OPENMS_RUNTIME_POST_EXCLUDE})
elseif(APPLE)
  set(OPENMS_RUNTIME_EXCLUDE "/usr/lib" "/System/")
  set(OPENMS_RUNTIME_POST_EXCLUDE "")
  # Backward compatibility
  set(EXCLUDE ${OPENMS_RUNTIME_EXCLUDE})
  set(POST_EXCLUDE ${OPENMS_RUNTIME_POST_EXCLUDE})
else()
  set(OPENMS_RUNTIME_EXCLUDE ".*/ld-linux-.*" ".*/linux-vdso.*" ".*/libm\\..*" ".*/libc\\..*" ".*/libpthread\\..*" ".*/libdl\\..*")
  set(OPENMS_RUNTIME_POST_EXCLUDE "")
  # Backward compatibility
  set(EXCLUDE ${OPENMS_RUNTIME_EXCLUDE})
  set(POST_EXCLUDE ${OPENMS_RUNTIME_POST_EXCLUDE})
endif()

# TODO check if we can reduce the permissions
install(RUNTIME_DEPENDENCY_SET OPENMS_DEPS
        DESTINATION ${INSTALL_LIB_DIR}
        PERMISSIONS
          OWNER_READ OWNER_WRITE OWNER_EXECUTE
          GROUP_READ GROUP_WRITE GROUP_EXECUTE
          WORLD_READ WORLD_WRITE WORLD_EXECUTE
        COMPONENT Dependencies
        PRE_EXCLUDE_REGEXES ${OPENMS_RUNTIME_EXCLUDE}
        POST_EXCLUDE_REGEXES ${OPENMS_RUNTIME_POST_EXCLUDE}
        DIRECTORIES $<TARGET_FILE_DIR:OpenMS>)

# Install GUI dependencies with same exclusion rules
if(WITH_GUI)
  install(RUNTIME_DEPENDENCY_SET OPENMS_GUI_DEPS
          DESTINATION ${INSTALL_LIB_DIR}
          PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_WRITE GROUP_EXECUTE
            WORLD_READ WORLD_WRITE WORLD_EXECUTE
          COMPONENT Dependencies
          PRE_EXCLUDE_REGEXES ${OPENMS_RUNTIME_EXCLUDE}
          POST_EXCLUDE_REGEXES ${OPENMS_RUNTIME_POST_EXCLUDE})
endif()

#install(RUNTIME_DEPENDENCY_SET TOPPView_DEPS) # I think without giving DESTINATION and COMPONENT it will be inferred
#install(RUNTIME_DEPENDENCY_SET TOPPAS_DEPS)
#...

########################################################### SEARCHENGINES
set(THIRDPARTY_COMPONENT_GROUP)
## populates the THIRDPARTY_COMPONENT_GROUP list
if(EXISTS ${SEARCH_ENGINES_DIRECTORY})
  ## TODO we could think about just recursing over subfolders
  install_thirdparty_folder("pwiz-bin")
  install_thirdparty_folder("Comet")
  install_thirdparty_folder("MSGFPlus")
  install_thirdparty_folder("XTandem")
  install_thirdparty_folder("LuciPHOr2")
  install_thirdparty_folder("SpectraST")
  install_thirdparty_folder("Sirius")
  install_thirdparty_folder("Percolator")
  install_thirdparty_folder("MaRaCluster")
  install_thirdparty_folder("ThermoRawFileParser")
endif()

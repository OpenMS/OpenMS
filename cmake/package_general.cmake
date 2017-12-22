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


# --------------------------------------------------------------------------
# general definitions used for building OpenMS packages
set(CPACK_PACKAGE_NAME "OpenMS")
set(CPACK_PACKAGE_VENDOR "OpenMS.de")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMS - A framework for mass spectrometry")
set(CPACK_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION_MAJOR "${OPENMS_PACKAGE_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${OPENMS_PACKAGE_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${OPENMS_PACKAGE_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMS-${CPACK_PACKAGE_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageDescriptionFile.cmake)
set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/License.txt)
set(CPACK_RESOURCE_FILE_WELCOME ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceWelcomeFile.txt)
set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceReadme.txt)

########################################################### Fixing dynamic dependencies
# Done on Windows via copying external and internal dlls to the install/bin/ folder
# Done on Mac via fixup_bundle for the GUI apps (TOPPView, TOPPAS) and via fix_mac_dependencies for the TOPP tools
# which recursively gathers dylds, copies them to install/lib/ and sets the install_name of the binaries to @executable_path/../lib
# Not done on Linux. Either install systemwide (omit CMAKE_INSTALL_PREFIX or set it to /usr/) or install and add the
# install/lib/ folder to the LD_LIBRARY_PATH

#install(CODE "
#  include(BundleUtilities)
#  GET_BUNDLE_ALL_EXECUTABLES(\${CMAKE_INSTALL_PREFIX}/${INSTALL_BIN_DIR} EXECS)
#  fixup_bundle(\"${EXECS}\" \"\" \"\${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}\")
#  " COMPONENT applications)

########################################################### SEARCHENGINES
set(THIRDPARTY_COMPONENT_GROUP)
## populates the THIRDPARTY_COMPONENT_GROUP list
if(EXISTS ${SEARCH_ENGINES_DIRECTORY})
  ## TODO we could think about just recursing over subfolders
  install_thirdparty_folder("Comet")
  install_thirdparty_folder("Fido")
  install_thirdparty_folder("MSGFPlus")
  install_thirdparty_folder("OMSSA")
  install_thirdparty_folder("XTandem")
  install_thirdparty_folder("LuciPHOr2")
  install_thirdparty_folder("MyriMatch")
  install_thirdparty_folder("SpectraST")
  install_thirdparty_folder("Sirius")
  install_thirdparty_folder("Percolator")
endif()

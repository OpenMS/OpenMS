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
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

include(CMakePackageConfigHelpers)

# a collection of wrapper for export functions that allows easier usage
# througout the OpenMS build system
set(_OPENMS_EXPORT_FILE "OpenMSTargets.cmake")

# clear list before we refill it
set(_OPENMS_EXPORT_TARGETS "" CACHE INTERNAL "List of targets that will be exported.")

macro(openms_register_export_target target_name)
  set(_OPENMS_EXPORT_TARGETS ${_OPENMS_EXPORT_TARGETS} ${target_name}
    CACHE INTERNAL "List of targets that will be exported.")
endmacro()

macro(openms_export_targets )

  # configure OpenMSConfig.cmake
  configure_package_config_file(
    "${OPENMS_HOST_DIRECTORY}/cmake/OpenMSConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/OpenMSConfig.cmake"
    INSTALL_DESTINATION ${INSTALL_LIB_DIR}/cmake/OpenMS
    PATH_VARS INSTALL_SHARE_DIR INSTALL_LIB_DIR INSTALL_DOC_DIR INSTALL_BIN_DIR
  )

  # write OpenMSConfigVersion.cmake
  write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/OpenMSConfigVersion.cmake"
    VERSION ${OPENMS_PACKAGE_VERSION}
    COMPATIBILITY SameMinorVersion
  )

  # create corresponding target file
  export(TARGETS ${_OPENMS_EXPORT_TARGETS}
         FILE ${OPENMS_HOST_BINARY_DIRECTORY}/${_OPENMS_EXPORT_FILE})

  # install the generated config file
  install_file(${PROJECT_BINARY_DIR}/OpenMSConfig.cmake
               ${INSTALL_CMAKE_DIR}
               cmake)

  # .. and ConfigVersion.cmake
  install_file(${PROJECT_BINARY_DIR}/OpenMSConfigVersion.cmake
               ${INSTALL_CMAKE_DIR}
               cmake)

  # register the package
  export(PACKAGE OpenMS)
endmacro()

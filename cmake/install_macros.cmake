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

# a collection of wrapper for install functions that allows easier usage
# througout the OpenMS build system

set(OPENMS_EXPORT_SET "OpenMSTargets")

#------------------------------------------------------------------------------
# Installs the library lib_target_name and all its headers set via
# set_target_properties(lib_target_name PROPERTIES PUBLIC_HEADER ${headers})
#
# @param lib_target_name The target name of the library that should be installed
macro(install_library lib_target_name)
  if ("${PACKAGE_TYPE}" STREQUAL "none")
    install(TARGETS ${lib_target_name}
      EXPORT ${OPENMS_EXPORT_SET}
      LIBRARY DESTINATION ${INSTALL_LIB_DIR}
      ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
      RUNTIME DESTINATION ${INSTALL_LIB_DIR}
      COMPONENT library)
  endif()
endmacro()

#------------------------------------------------------------------------------
# Installs the given headers.
#
# @param header_list List of headers to install
macro(install_headers header_list component)
  foreach(_header ${header_list})
    set(_relative_header_path)

    get_filename_component(_target_path ${_header} PATH)
    if ("${_target_path}" MATCHES "^${PROJECT_BINARY_DIR}.*")
      # is generated bin header
      string(REPLACE "${PROJECT_BINARY_DIR}/include/OpenMS" "" _relative_header_path "${_target_path}")
    else()
      # is source header -> strip include/OpenMS
      string(REPLACE "include/OpenMS" "" _relative_header_path "${_target_path}")
    endif()

    # install the header
    install(FILES ${_header}
            # note the missing slash, we need this for file directly located in
            # include/OpenMS (e.g., config.h)
            DESTINATION ${INSTALL_INCLUDE_DIR}/OpenMS${_relative_header_path}
            COMPONENT ${component}_headers)
  endforeach()
endmacro()

#------------------------------------------------------------------------------
# Installs the tool tool_target_name
# @param tool_target_name The target name of the tool that should be installed
macro(install_tool tool_target_name)
  if ("${PACKAGE_TYPE}" STREQUAL "none")
    install(TARGETS ${tool_target_name}
      RUNTIME DESTINATION ${INSTALL_BIN_DIR}
      BUNDLE DESTINATION ${INSTALL_BIN_DIR}
      COMPONENT applications)
  endif()
endmacro()

#------------------------------------------------------------------------------
# Installs a given directory
# @param directory The directory to install
# @param destination The destination (relative to the prefix) where it should be installed
# @param component The component to which to the directory belongs
macro(install_directory directory destination component)
  if ("${PACKAGE_TYPE}" STREQUAL "none")
    install(DIRECTORY ${directory}
      DESTINATION ${destination}
      COMPONENT ${component})
  endif()
endmacro()

#------------------------------------------------------------------------------
# Installs a given file
# @param directory The file to install
# @param destination The destination (relative to the prefix) where it should be installed
# @param component The component to which to the file belongs
macro(install_file file destination component)
  if ("${PACKAGE_TYPE}" STREQUAL "none")
    install(FILES ${file}
      DESTINATION ${destination}
      COMPONENT ${component})
  endif()
endmacro()

#------------------------------------------------------------------------------
# Execute the given code while executing the install target
# @param code_snippet The code to execute
# @param component The component to which the code will be associated
macro(install_code code_snippet component)
  if ("${PACKAGE_TYPE}" STREQUAL "none")
    install(CODE ${code_snippet}
            COMPONENT ${component})
  endif()
endmacro()

#------------------------------------------------------------------------------
# Installs the exported target information
macro(install_export_targets )
  if ("${PACKAGE_TYPE}" STREQUAL "none")
    install(EXPORT ${OPENMS_EXPORT_SET}
            DESTINATION ${INSTALL_SHARE_DIR}/cmake
            COMPONENT share)
  endif()
endmacro()

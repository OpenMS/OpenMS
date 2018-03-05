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

# a collection of wrapper for install functions that allows easier usage
# throughout the OpenMS build system

set(OPENMS_EXPORT_SET "OpenMSTargets")

#------------------------------------------------------------------------------
# Installs the library lib_target_name and all its headers set via
# set_target_properties(lib_target_name PROPERTIES PUBLIC_HEADER ${headers})
#
# @param lib_target_name The target name of the library that should be installed
macro(install_library lib_target_name)
    install(TARGETS ${lib_target_name}
      EXPORT ${OPENMS_EXPORT_SET}
      LIBRARY DESTINATION ${INSTALL_LIB_DIR} COMPONENT library
      ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT library
      RUNTIME DESTINATION ${INSTALL_LIB_DIR} COMPONENT library
      )
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
    install(TARGETS ${tool_target_name}
      RUNTIME DESTINATION ${INSTALL_BIN_DIR} COMPONENT Applications
      BUNDLE DESTINATION ${INSTALL_BIN_DIR} COMPONENT Applications
      )
endmacro()

#------------------------------------------------------------------------------
# Installs a given directory
# @param directory The directory to install
# @param destination The destination (relative to the prefix) where it should be installed
# @param component The component to which to the directory belongs
macro(install_directory directory destination component)
    install(DIRECTORY ${directory}
      DESTINATION ${destination}
      COMPONENT ${component}
      FILE_PERMISSIONS      OWNER_WRITE OWNER_READ
                            GROUP_READ
                            WORLD_READ
      DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                            GROUP_EXECUTE GROUP_READ
                            WORLD_EXECUTE WORLD_READ
      REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
      REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
        )
endmacro()

#------------------------------------------------------------------------------
# Installs a given file
# @param directory The file to install
# @param destination The destination (relative to the prefix) where it should be installed
# @param component The component to which to the file belongs
macro(install_file file destination component)
    install(FILES ${file}
      DESTINATION ${destination}
      COMPONENT ${component})
endmacro()

#------------------------------------------------------------------------------
# Execute the given code while executing the install target
# @param code_snippet The code to execute
# @param component The component to which the code will be associated
macro(install_code code_snippet component)
    install(CODE ${code_snippet}
            COMPONENT ${component})
endmacro()

#------------------------------------------------------------------------------
# Installs the exported target information
macro(install_export_targets )
    install(EXPORT ${OPENMS_EXPORT_SET}
            DESTINATION ${INSTALL_SHARE_DIR}/cmake
            COMPONENT share)
endmacro()

#------------------------------------------------------------------------------
# Installs Thirdparty folders with executables
macro(install_thirdparty_folder foldername)
  if(EXISTS ${SEARCH_ENGINES_DIRECTORY}/${foldername})
    install(DIRECTORY             ${SEARCH_ENGINES_DIRECTORY}/${foldername}
            DESTINATION           ${INSTALL_SHARE_DIR}/THIRDPARTY
            COMPONENT             ${foldername}
            FILE_PERMISSIONS      OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                                  GROUP_READ GROUP_EXECUTE
                                  WORLD_READ WORLD_EXECUTE
            REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
            REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
            )
    list(APPEND THIRDPARTY_COMPONENT_GROUP ${foldername})
  endif()
endmacro()

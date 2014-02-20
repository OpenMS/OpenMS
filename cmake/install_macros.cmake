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

# a collection of wrapper for install functions that allows easier usage
# througout the OpenMS build system

#------------------------------------------------------------------------------
# customize output folders according to platform
if (WIN32)
	set(OPENMS_LIB_INSTALL_PATH "bin")
else()	## Linux & MacOS
  set(OPENMS_LIB_INSTALL_PATH "lib")
endif()

#------------------------------------------------------------------------------
# Installs the library lib_target_name
# @param lib_target_name The target name of the library that should be installed
macro(install_library lib_target_name)
  if ("${PACKAGE_TYPE}" STREQUAL "none")
    install(TARGETS ${lib_target_name}
      LIBRARY DESTINATION ${OPENMS_LIB_INSTALL_PATH}
      ARCHIVE DESTINATION ${OPENMS_LIB_INSTALL_PATH}
    	RUNTIME DESTINATION ${OPENMS_LIB_INSTALL_PATH}
      COMPONENT library)
  endif()
endmacro()

#------------------------------------------------------------------------------
# Installs the tool tool_target_name
# @param tool_target_name The target name of the tool that should be installed
macro(install_tool tool_target_name)
  if ("${PACKAGE_TYPE}" STREQUAL "none")
    install(TARGETS ${tool_target_name}
      RUNTIME DESTINATION bin
      BUNDLE DESTINATION bin
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
# Executes the given code on execution of install
# @param directory The file to install
macro(execute_on_install code)
  if ("${PACKAGE_TYPE}" STREQUAL "none")
    install(CODE "${code}"
      COMPONENT 01) # we give this weird component name to ensure that it is executed before anything else
  endif()
endmacro()
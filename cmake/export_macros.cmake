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

# a collection of wrapper for export functions that allows easier usage
# througout the OpenMS build system
set(_OPENMS_EXPORT_FILE "${OPENMS_HOST_BINARY_DIRECTORY}/OpenMSTargets.cmake")

# clear list before we refill it
set(_OPENMS_EXPORT_TARGETS "" CACHE INTERNAL "List of targets that will be exported.")

macro(openms_register_export_target target_name)
  set(_OPENMS_EXPORT_TARGETS ${_OPENMS_EXPORT_TARGETS} ${target_name}
    CACHE INTERNAL "List of targets that will be exported.")

	get_target_property(_TARGET_COMPILE_DEFINITIONS
                ${target_name}
                COMPILE_DEFINITIONS)

  message(STATUS "${target_name} cmp_defs: ${_TARGET_COMPILE_DEFINITIONS}")
endmacro()

macro(openms_export_targes )
  set(_EXPORT_INCLUDE_BLOCK "")

  # we also need to export the corresponding include directories
  foreach(_target ${_OPENMS_EXPORT_TARGETS})
    # check if we have a corresponding include_dir variable
    if(NOT DEFINED ${_target}_INCLUDE_DIRECTORIES)
      message(FATAL_ERROR "Please provide the matching include directory variable ${_target}_INCLUDE_DIRECTORIES for export target ${_target}")
    endif()

    # extend include block
    set(_EXPORT_INCLUDE_BLOCK "set(${_target}_INCLUDE_DIRECTORIES ${${_target}_INCLUDE_DIRECTORIES})\n${_EXPORT_INCLUDE_BLOCK}")
  endforeach()

  # configure OpenMSConfig.cmake
	configure_file(
		"${OPENMS_HOST_DIRECTORY}/cmake/OpenMSConfig.cmake.in"
		"${PROJECT_BINARY_DIR}/OpenMSConfig.cmake"
		@ONLY
	)

  # create corresponding target file
  export(TARGETS ${_OPENMS_EXPORT_TARGETS} FILE ${_OPENMS_EXPORT_FILE})

  # register the package
	export(PACKAGE OpenMS)
endmacro()

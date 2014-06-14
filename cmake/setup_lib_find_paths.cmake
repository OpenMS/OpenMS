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
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# This cmake file only handles the customization of internal path variables
# where the build system expects to find external libraries. Note that the
# actual libraries are found in the CMake files of the individual componenents.
#------------------------------------------------------------------------------

set(CONTRIB_CUSTOM_DIR CACHE DOC "DEPRECATED: Please use CMAKE_FIND_ROOT_PATH instead! User defined location of contrib dir. If left empty we assume the contrib to be in OpenMS/contrib!")
set(CONTRIB_DIR ${OPENMS_HOST_DIRECTORY}/contrib/ CACHE INTERNAL "Final contrib path after looking at CMAKE_FIND_ROOT_PATH. Defaults to OpenMS/contrib")

if("${CMAKE_FIND_ROOT_PATH}" STREQUAL "")
	if(NOT "${CONTRIB_CUSTOM_DIR}" STREQUAL "")
		message("Please do no longer use -DCONTRIB_CUSTOM_DIR. This option is deprecated. Please use -DCMAKE_FIND_ROOT_PATH instead.")
		list(INSERT CONTRIB_DIR 0 ${CONTRIB_CUSTOM_DIR})
  else()
		# Append a few unusual default search directories for convenience
		# if no FIND ROOT PATH was specified.
		list(APPEND CONTRIB_DIR /opt/local /usr/local /usr)
	endif()
endif()

#------------------------------------------------------------------------------
# Append all contrib dirs to the (potentially empty) FIND_ROOT_PATH.
# This will be the final search order used by regular CMAKE modules (by default)
# and by the OpenMS macros (via CONTRIB_DIR).
list(APPEND CMAKE_FIND_ROOT_PATH "${CONTRIB_DIR}")
set(TMP "")
foreach(CUSTOM_PATH ${CMAKE_FIND_ROOT_PATH})
  get_filename_component(ABS_PATH ${CUSTOM_PATH} ABSOLUTE)
	list(APPEND TMP ${ABS_PATH})
endforeach()
set(CMAKE_FIND_ROOT_PATH "${TMP}")
set(CONTRIB_DIR "${CMAKE_FIND_ROOT_PATH}")

set(CONTRIB_INCLUDE_DIR "" CACHE INTERNAL "contrib include dir")
set(CONTRIB_LIB_DIR "" CACHE INTERNAL "contrib lib dir")
foreach(CONTRIB_PATH ${CONTRIB_DIR})
  list(APPEND CONTRIB_INCLUDE_DIR "${CONTRIB_PATH}/include/")
  list(APPEND CONTRIB_LIB_DIR "${CONTRIB_PATH}/lib/")
endforeach()
message(STATUS "===========================================================================")
message(STATUS "CMake find root path: ${CMAKE_FIND_ROOT_PATH}")
message(STATUS "Contrib search directories:  ${CONTRIB_DIR}")
message(STATUS "Contrib library directories: ${CONTRIB_LIB_DIR}")
message(STATUS "Contrib include directories: ${CONTRIB_INCLUDE_DIR}")
message(STATUS "===========================================================================")

#------------------------------------------------------------------------------
# Ensure Qt includes it's libs as SYSTEM
set(QT_INCLUDE_DIRS_NO_SYSTEM Off)

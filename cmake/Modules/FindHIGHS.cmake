# - Try to find HIGHS
# Once done this will define
#
#  HIGHS_FOUND - system has HIGHS
#  HIGHS_INCLUDE_DIRS - the HIGHS include directory
#  HIGHS_LIBRARIES - Link these to use HIGHS
#
# Based on FindGLPK.cmake and adapted for HIGHS
#
#=============================================================================
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
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN if
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# required scripts
include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)

# hint from the user
set(HIGHS_ROOT_DIR "" CACHE PATH "HIGHS root directory")

find_path(HIGHS_INCLUDE_DIR
  Highs.h
  PATHS ${HIGHS_ROOT_DIR}/include
  PATH_SUFFIXES highs
)

if (NOT HIGHS_LIBRARIES)
	# find library and (optionally) the corresponding debug version
	find_library(HIGHS_LIBRARY_RELEASE
		NAMES highs
		PATHS ${HIGHS_ROOT_DIR}/lib
	)

	find_library(HIGHS_LIBRARY_DEBUG
		NAMES highsd
		PATHS ${HIGHS_ROOT_DIR}/lib
	)

	SELECT_LIBRARY_CONFIGURATIONS(HIGHS)
endif()

if(HIGHS_INCLUDE_DIR AND HIGHS_LIBRARY)
  # Try to extract version from header if available
  if(EXISTS ${HIGHS_INCLUDE_DIR}/Highs.h)
    file(READ ${HIGHS_INCLUDE_DIR}/Highs.h HIGHS_HIGHS_H)
    
    string(REGEX MATCH "define[ ]+HIGHS_VERSION_MAJOR[ ]+[0-9]+" HIGHS_MAJOR_VERSION_LINE "${HIGHS_HIGHS_H}")
    string(REGEX REPLACE "define[ ]+HIGHS_VERSION_MAJOR[ ]+([0-9]+)" "\\1" HIGHS_VERSION_MAJOR "${HIGHS_MAJOR_VERSION_LINE}")
    
    string(REGEX MATCH "define[ ]+HIGHS_VERSION_MINOR[ ]+[0-9]+" HIGHS_MINOR_VERSION_LINE "${HIGHS_HIGHS_H}")
    string(REGEX REPLACE "define[ ]+HIGHS_VERSION_MINOR[ ]+([0-9]+)" "\\1" HIGHS_VERSION_MINOR "${HIGHS_MINOR_VERSION_LINE}")
    
    string(REGEX MATCH "define[ ]+HIGHS_VERSION_PATCH[ ]+[0-9]+" HIGHS_PATCH_VERSION_LINE "${HIGHS_HIGHS_H}")
    string(REGEX REPLACE "define[ ]+HIGHS_VERSION_PATCH[ ]+([0-9]+)" "\\1" HIGHS_VERSION_PATCH "${HIGHS_PATCH_VERSION_LINE}")
    
    if(HIGHS_VERSION_MAJOR AND HIGHS_VERSION_MINOR)
      set(HIGHS_VERSION_STRING "${HIGHS_VERSION_MAJOR}.${HIGHS_VERSION_MINOR}")
      if(HIGHS_VERSION_PATCH)
        set(HIGHS_VERSION_STRING "${HIGHS_VERSION_STRING}.${HIGHS_VERSION_PATCH}")
      endif()
    endif()
  endif()
  
  # Set to true if we found the library
  set(HIGHS_PROPER_VERSION_FOUND TRUE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HIGHS DEFAULT_MSG HIGHS_LIBRARY HIGHS_INCLUDE_DIR HIGHS_PROPER_VERSION_FOUND)

if(HIGHS_FOUND)
  set(HIGHS_INCLUDE_DIRS ${HIGHS_INCLUDE_DIR})
  set(HIGHS_LIBRARIES ${HIGHS_LIBRARY})
  if(NOT TARGET HIGHS::HIGHS)
    add_library(HIGHS::HIGHS UNKNOWN IMPORTED)
    set_property(TARGET HIGHS::HIGHS PROPERTY IMPORTED_LOCATION "${HIGHS_LIBRARY_RELEASE}")
    if(HIGHS_LIBRARY_DEBUG)
      set_property(TARGET HIGHS::HIGHS PROPERTY IMPORTED_LOCATION_DEBUG "${HIGHS_LIBRARY_DEBUG}")
    endif()
    set_property(TARGET HIGHS::HIGHS PROPERTY INCLUDE_DIRECTORIES "${HIGHS_INCLUDE_DIR}")
    set_property(TARGET HIGHS::HIGHS PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${HIGHS_INCLUDE_DIR}")
  endif()
endif(HIGHS_FOUND)

mark_as_advanced(HIGHS_LIBRARY HIGHS_INCLUDE_DIR)

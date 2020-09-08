# - Try to find Xerces-C++
# Once done this will define
#
#  XercesC_FOUND - system has XercesC
#  XercesC_INCLUDE_DIR - the XercesC include directory
#  XercesC_LIBRARIES - Link these to use XercesC
#  XercesC_VERSION_STRING - the version of XercesC found
#
# Inspired by Ben Morgan, <Ben.Morgan@warwick.ac.uk>
# http://geant4.cern.ch/support/source/geant4/cmake/Modules/FindXercesC.cmake
#
#=============================================================================
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
#=============================================================================

# additional search paths
set(_XercesC_PATHS
  "[HKEY_CURRENT_USER\\software\\xerces-c\\src]"
  "[HKEY_CURRENT_USER\\xerces-c\\src]"
)

set(_XercesC_INCLUDE_TARGET "xercesc/util/XercesVersion.hpp")

# Find Xerce-C include path
find_path(
    XercesC_INCLUDE_DIRS
    PATHS ${_XercesC_PATHS}
    NAMES ${_XercesC_INCLUDE_TARGET}
)

# Find the xerces libraries
if (NOT XercesC_LIBRARIES)
    ## The NAMES_PER_DIR option will make sure that the PATHS are the "outer for loop" when searching for the libraries.
	## We want that because we put the contrib as the first search path usually.
    find_library(XercesC_LIBRARY_RELEASE NAMES xerces-c xerces-c_3 xerces-c_3_1 xerces-c-3.1 xerces-c_3_2 xerces-c-3.2 NAMES_PER_DIR ${_XercesC_PATHS} PATH_SUFFIXES lib)
    find_library(XercesC_LIBRARY_DEBUG NAMES xerces-c_3D xerces-c_3_1D xerces-c-3.1D xerces-c_3_2D xerces-c-3.2D NAMES_PER_DIR ${_XercesC_PATHS} PATH_SUFFIXES lib)

    include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)
    select_library_configurations(XercesC)
endif ()

# identify xerces version
if (XercesC_INCLUDE_DIRS AND EXISTS "${XercesC_INCLUDE_DIRS}/${_XercesC_INCLUDE_TARGET}")
  file(STRINGS "${XercesC_INCLUDE_DIRS}/${_XercesC_INCLUDE_TARGET}" _XercesC_H REGEX "^#define XERCES_VERSION_.* [0-9]+")
  #define XERCES_VERSION_MAJOR 3
  string(REGEX REPLACE ".*\#define XERCES_VERSION_MAJOR ([0-9]+).*" "\\1" XercesC_VERSION_MAJOR "${_XercesC_H}")
  #define XERCES_VERSION_MINOR 1
  string(REGEX REPLACE ".*\#define XERCES_VERSION_MINOR ([0-9]+).*" "\\1" XercesC_VERSION_MINOR "${_XercesC_H}")
  #define XERCES_VERSION_REVISION 1
  string(REGEX REPLACE ".*\#define XERCES_VERSION_REVISION ([0-9]+).*" "\\1" XercesC_VERSION_REVISION "${_XercesC_H}")

  set(XercesC_VERSION_STRING "${XercesC_VERSION_MAJOR}.${XercesC_VERSION_MINOR}.${XercesC_VERSION_REVISION}")
endif ()

# handle the QUIETLY and REQUIRED arguments and set XercesC_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XercesC
                                  REQUIRED_VARS XercesC_LIBRARIES XercesC_INCLUDE_DIRS
                                  VERSION_VAR XercesC_VERSION_STRING)

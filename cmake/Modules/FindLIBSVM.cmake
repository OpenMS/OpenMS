# - Try to find LIBSVM
# Once done this will define
#
#  LIBSVM_FOUND - system has LIBSVM
#  LIBSVM_INCLUDE_DIR - the LIBSVM include directory
#  LIBSVM_LIBRARIES - Link these to use LIBSVM
#  LIBSVM_VERSION_STRING - the version of LIBSVM found
#
# Inspired by Julien Schueller <schueller at phimeca dot com>
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
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# set LIBSVM_INCLUDE_DIR
find_path (LIBSVM_INCLUDE_DIR NAMES svm.h PATH_SUFFIXES libsvm libsvm-3.1/libsvm DOC "LibSVM include directory" )

# set LIBSVM_INCLUDE_DIRS
if (NOT LIBSVM_INCLUDE_DIRS)
  set (LIBSVM_INCLUDE_DIRS ${LIBSVM_INCLUDE_DIR})
endif ()

# extract version
set (LIBSVM_MAJOR_VERSION 0)
set (LIBSVM_MINOR_VERSION 0)
set (LIBSVM_SUBMINOR_VERSION 0)
if (LIBSVM_INCLUDE_DIR)
  # LIBSVM_VERSION macro defined in svm.h since version 2.8.9
  file (STRINGS "${LIBSVM_INCLUDE_DIR}/svm.h" _VERSION_STRING REGEX ".*LIBSVM_VERSION.*")
  if (_VERSION_STRING)
    string (REGEX REPLACE ".*_VERSION[ ]+([0-9]+)" "\\1" _VERSION_NUMBER "${_VERSION_STRING}")
    math (EXPR LIBSVM_MAJOR_VERSION "${_VERSION_NUMBER} / 100")
    math (EXPR LIBSVM_MINOR_VERSION "(${_VERSION_NUMBER} % 100 ) / 10")
    math (EXPR LIBSVM_SUBMINOR_VERSION "${_VERSION_NUMBER} % 10")
  endif ()
endif ()
set (LIBSVM_VERSION "${LIBSVM_MAJOR_VERSION}.${LIBSVM_MINOR_VERSION}.${LIBSVM_SUBMINOR_VERSION}")

# find LIBSVM_LIBRARY
find_library (LIBSVM_LIBRARY_RELEASE NAMES svm libsvm DOC "LibSVM library location" )
find_library (LIBSVM_LIBRARY_DEBUG NAMES svmd libsvmd DOC "LibSVM library location" )

include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)
select_library_configurations(LIBSVM)

if(NOT TARGET LibSVM::LibSVM)
  add_library(LibSVM::LibSVM UNKNOWN IMPORTED) # TODO we could try to infer shared/static
  set_property(TARGET LibSVM::LibSVM PROPERTY IMPORTED_LOCATION "${LIBSVM_LIBRARY_RELEASE}")
  set_property(TARGET LibSVM::LibSVM PROPERTY IMPORTED_LOCATION_DEBUG "${LIBSVM_LIBRARY_DEBUG}")
  set_property(TARGET LibSVM::LibSVM PROPERTY INCLUDE_DIRECTORIES "${LIBSVM_INCLUDE_DIR}")
  set_property(TARGET LibSVM::LibSVM PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${LIBSVM_INCLUDE_DIR}")
endif()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (LIBSVM
                                  REQUIRED_VARS LIBSVM_LIBRARIES LIBSVM_INCLUDE_DIRS
                                  VERSION_VAR LIBSVM_VERSION)
mark_as_advanced (
  LIBSVM_LIBRARIES
  LIBSVM_INCLUDE_DIR
  LIBSVM_INCLUDE_DIRS
  LIBSVM_VERSION
)

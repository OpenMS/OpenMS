# - Try to find SQLITE
# Once done this will define
#
#  SQLITE_FOUND - system has SQLITE
#  SQLITE_INCLUDE_DIR - the SQLITE include directory
#  SQLITE_LIBRARIES - Link these to use SQLITE
#  SQLITE_VERSION_STRING - the version of SQLITE found
#
# Inspired by Julien Schueller <schueller at phimeca dot com>
#
#=============================================================================
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
#=============================================================================

# set SQLITE_INCLUDE_DIR
find_path(SQLITE_INCLUDE_DIR sqlite3.h PATH_SUFFIXES ./include/sqlite)

# extract version
if (SQLITE_INCLUDE_DIR)
  file (STRINGS "${SQLITE_INCLUDE_DIR}/sqlite3.h" _VERSION_STRING REGEX ".*SQLITE_VERSION.*")
endif()

# find SQLITE_LIBRARY
find_library (SQLITE_LIBRARY NAMES sqlite3 SQLITE "SQLITE library location" )

include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)
select_library_configurations(SQLITE)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SQLITE
                                  REQUIRED_VARS SQLITE_LIBRARY SQLITE_INCLUDE_DIR
                                  VERSION_VAR SQLITE_VERSION)
mark_as_advanced (
  SQLITE_LIBRARY
  SQLITE_INCLUDE_DIR
  SQLITE_VERSION
)

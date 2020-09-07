# - Try to find COIN
# Once done this will define
#
#  COIN_FOUND - system has Coin
#  COIN_INCLUDE_DIRS - the Coin include directory
#  COIN_LIBRARIES - Link these to use Coin
#
# Inspired by LEMON's FindCOIN
# https://lemon.cs.elte.hu/trac/lemon/browser/lemon/cmake/FindCOIN.cmake
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
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN if
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

# required scripts
include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)

# hint from the user
set(COIN_ROOT_DIR "" CACHE PATH "COIN root directory")

# find the coin include directory
find_path(COIN_INCLUDE_DIR coin/CoinUtilsConfig.h
  HINTS ${COIN_ROOT_DIR}/include
)

# helper macro to find specific coind sub-libraries
macro(_coin_find_lib _libname _lib_file_names _lib_file_names_debug)
  if(NOT COIN_${_libname})
    # find release version
    find_library(COIN_${_libname}_LIBRARY_RELEASE
      NAMES ${_lib_file_names}
      HINTS ${COIN_ROOT_DIR}/lib/coin
      HINTS ${COIN_ROOT_DIR}/lib
    )
    # .. and debug version
    find_library(COIN_${_libname}_LIBRARY_DEBUG
      NAMES ${_lib_file_names_debug}
      HINTS ${COIN_ROOT_DIR}/lib/coin
      HINTS ${COIN_ROOT_DIR}/lib
    )

    # create final library to be exported
    select_library_configurations(COIN_${_libname})
  endif()
endmacro()

# actually find the required libs
_coin_find_lib("CBC" "libCbc;Cbc" "libCbcd;Cbc")
_coin_find_lib("CGL" "libCgl;Cgl" "libCgld;Cgl")
_coin_find_lib("CLP" "libClp;Clp" "libClpd;Clp")
_coin_find_lib("COIN_UTILS" "libCoinUtils;CoinUtils" "libCoinUtilsd;CoinUtils")
_coin_find_lib("OSI" "libOsi;Osi" "libOsid;Osi")
_coin_find_lib("OSI_CLP" "libOsiClp;OsiClp" "libOsiClpd;OsiClp")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(COIN DEFAULT_MSG
  COIN_INCLUDE_DIR
  COIN_CBC_LIBRARY
  COIN_CGL_LIBRARY
  COIN_CLP_LIBRARY
  COIN_COIN_UTILS_LIBRARY
  COIN_OSI_LIBRARY
  COIN_OSI_CLP_LIBRARY
)

# export the libraries
if(COIN_FOUND)
  set(COIN_INCLUDE_DIRS ${COIN_INCLUDE_DIR})
  set(COIN_LIBRARIES 
    ${COIN_CBC_LIBRARY}
    ${COIN_CGL_LIBRARY}
    ${COIN_CLP_LIBRARY}
    ${COIN_COIN_UTILS_LIBRARY}
    ${COIN_OSI_LIBRARY}
    ${COIN_OSI_CLP_LIBRARY}
  )
endif(COIN_FOUND)

mark_as_advanced(
  COIN_INCLUDE_DIR
  COIN_CBC_LIBRARY
  COIN_CGL_LIBRARY
  COIN_CLP_LIBRARY
  COIN_COIN_UTILS_LIBRARY
  COIN_OSI_LIBRARY
  COIN_OSI_CLP_LIBRARY
)

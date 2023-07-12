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
# ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

# find for vcpkg
find_path(COIN_INCLUDE_DIR coin-or/CoinUtilsConfig.h
        HINTS
        ${COIN_ROOT_DIR}/include
        )

if (NOT COIN_INCLUDE_DIR)
  # find the coin include directory from contrib or system
  find_path(COIN_INCLUDE_DIR coin-or/CoinUtilsConfig.h coin/CoinUtilsConfig.h coinutils/coin/CoinUtilsConfig.h
          HINTS
          ${COIN_ROOT_DIR}/include
          )
  if (COIN_INCLUDE_DIR)
    set(CF_COIN_INCLUDE_SUBDIR_DEF 1 CACHE BOOL "If the subdir for including coin-or headers is coin (1) or coin-or (undefined).")
  endif()
endif()

# helper macro to find specific coin sub-libraries
macro(_coin_find_lib _libname _lib_file_names _lib_file_names_debug)
  if(NOT COIN_${_libname})
    string(TOLOWER ${_libname} _libnamelower)
    if(_libnamelower STREQUAL "osi_clp")
      set(_libnamelower "clp")
    endif()

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

    if(EXISTS "${COIN_INCLUDE_DIR}/${_libnamelower}/coin")
      ## Unfortunately coin uses a weird include model and includes headers from different subpackages
      ## via include<File.hpp> instead of include<coin/File.hpp>
      set(_libspecific_include "${COIN_INCLUDE_DIR}/${_libnamelower}" "${COIN_INCLUDE_DIR}/${_libnamelower}/coin")
    else()
      set(_libspecific_include)
    endif()

    # create final library to be exported
    select_library_configurations(COIN_${_libname})
    if(NOT TARGET COIN_${_libname})
      add_library(CoinOR::${_libname} UNKNOWN IMPORTED) # TODO we could try to infer shared/static
      set_property(TARGET CoinOR::${_libname} PROPERTY IMPORTED_LOCATION "${COIN_${_libname}_LIBRARY_RELEASE}")
      set_property(TARGET CoinOR::${_libname} PROPERTY IMPORTED_LOCATION_DEBUG "${COIN_${_libname}_LIBRARY_DEBUG}")
      set_property(TARGET CoinOR::${_libname} PROPERTY INCLUDE_DIRECTORIES "${COIN_INCLUDE_DIR}" ${_libspecific_include})
      set_property(TARGET CoinOR::${_libname} PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${COIN_INCLUDE_DIR}" ${_libspecific_include})
      # TODO there are probably dependencies across the single libs but we don't care for now. We always include all.
      set_property(TARGET CoinOR::CoinOR APPEND PROPERTY INTERFACE_LINK_LIBRARIES CoinOR::${_libname})
    endif()
  endif()
endmacro()


if(NOT TARGET CoinOR::CoinOR)
  add_library(CoinOR::CoinOR INTERFACE IMPORTED)
  if (VCPKG_TOOLCHAIN)
    # Currently coin-or from vcpkg requires BLAS and LAPACK
    find_package(BLAS)
    find_package(LAPACK)
    target_link_libraries(CoinOR::CoinOR INTERFACE BLAS LAPACK)
  endif()
endif()

# actually find the required libs and add them as dependencies
# to the parent CoinOR::CoinOR interface target
_coin_find_lib("CBC" "libCbc;Cbc" "libCbcd;Cbc")
_coin_find_lib("CGL" "libCgl;Cgl" "libCgld;Cgl")
_coin_find_lib("CLP" "libClp;Clp" "libClpd;Clp")
_coin_find_lib("COINUTILS" "libCoinUtils;CoinUtils" "libCoinUtilsd;CoinUtils")
_coin_find_lib("OSI" "libOsi;Osi" "libOsid;Osi")
_coin_find_lib("OSI_CLP" "libOsiClp;OsiClp" "libOsiClpd;OsiClp")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(COIN DEFAULT_MSG
  COIN_INCLUDE_DIR
  COIN_CBC_LIBRARY
  COIN_CGL_LIBRARY
  COIN_CLP_LIBRARY
  COIN_COINUTILS_LIBRARY
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
    ${COIN_COINUTILS_LIBRARY}
    ${COIN_OSI_LIBRARY}
    ${COIN_OSI_CLP_LIBRARY}
  )
endif(COIN_FOUND)

mark_as_advanced(
  COIN_INCLUDE_DIR
  COIN_CBC_LIBRARY
  COIN_CGL_LIBRARY
  COIN_CLP_LIBRARY
  COIN_COINUTILS_LIBRARY
  COIN_OSI_LIBRARY
  COIN_OSI_CLP_LIBRARY
)

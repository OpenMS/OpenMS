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
set(COIN_ROOT_DIR "" CACHE PATH "COIN root directory")

# find for vcpkg
find_path(COIN_VCPKG_INCLUDE_DIR coin-or/CoinUtilsConfig.h
  HINTS
  ${COIN_ROOT_DIR}/include
)

# find for contrib and system
find_path(COIN_SYS_INCLUDE_DIR coin/CoinUtilsConfig.h coinutils/coin/CoinUtilsConfig.h
  HINTS
  ${COIN_ROOT_DIR}/include
)

if (COIN_SYS_INCLUDE_DIR)
  set(COIN_INCLUDE_DIR ${COIN_SYS_INCLUDE_DIR})
  set(OPENMS_HAS_COIN_INCLUDE_SUBDIR_IS_COIN 1 CACHE BOOL "If the subdir for including coin-or headers is 'coin' (1) or 'coin-or' (undefined).")
elseif (COIN_VCPKG_INCLUDE_DIR)
  set(COIN_INCLUDE_DIR ${COIN_VCPKG_INCLUDE_DIR})
  unset(OPENMS_HAS_COIN_INCLUDE_SUBDIR_IS_COIN)
endif() # find_package_handle_standard_args will handle missingness

# helper macro to find specific coin sub-libraries
macro(_coin_find_lib _libname _libname_camel _lib_file_names _lib_file_names_debug)
  if(NOT COIN_${_libname})
    if (${_libname} STREQUAL "OSI_CLP")
      set(FOLDER "clp")
      set(HNAME "OsiClpSolverInterface.hpp")
    else()
      string(TOLOWER ${_libname_camel} FOLDER)
      set(HNAME ${_libname_camel}Config.h)
    endif()

    # find release version
    find_library(COIN_${_libname}_LIBRARY_RELEASE
      NAMES ${_lib_file_names}
      HINTS ${COIN_ROOT_DIR}/lib/coin
            ${COIN_ROOT_DIR}/lib
    )
    # .. and debug version
    find_library(COIN_${_libname}_LIBRARY_DEBUG
      NAMES ${_lib_file_names_debug}
      HINTS ${COIN_ROOT_DIR}/lib/coin
      HINTS ${COIN_ROOT_DIR}/lib
    )

    find_path(
      ${_libname}_INCLUDE_DIR
      NAMES 
        ${FOLDER}/coin/${HNAME} 
    )

    set (INCLUDE_DIRS)
    if (EXISTS "${${_libname}_INCLUDE_DIR}/${FOLDER}")
      # we (in OpenMS) include the coin headers with e.g., coin/ClpInterface.h
      list(APPEND INCLUDE_DIRS "${${_libname}_INCLUDE_DIR}/${FOLDER}")
      # while coin-or includes headers from other sublibraries with just e.g., ClpModel.h
      list(APPEND INCLUDE_DIRS "${${_libname}_INCLUDE_DIR}/${FOLDER}/coin")
    endif()

    if (EXISTS ${COIN_INCLUDE_DIR})
      list(APPEND INCLUDE_DIRS ${COIN_INCLUDE_DIR})
    endif()

    # Before, we handled the vcpkg specific sublibrary structure as follows. This should be handled equally
    # by the more generic searching above (which also works for brew). But needs to be tested. I don't use
    # vcpkg currently
    # _libnamelower = FOLDER

    #if(EXISTS "${COIN_INCLUDE_DIR}/${_libnamelower}/coin")
    #  set(_libspecific_include "${COIN_INCLUDE_DIR}/${_libnamelower}" "${COIN_INCLUDE_DIR}/${_libnamelower}/coin")
    #else()
    #  set(_libspecific_include)
    #endif()

    # create final library to be exported
    select_library_configurations(COIN_${_libname})
    if(NOT TARGET COIN_${_libname})
      add_library(CoinOR::${_libname} UNKNOWN IMPORTED) # TODO we could try to infer shared/static instead of UNKNOWN
      set_property(TARGET CoinOR::${_libname} PROPERTY IMPORTED_LOCATION "${COIN_${_libname}_LIBRARY_RELEASE}")
      set_property(TARGET CoinOR::${_libname} PROPERTY IMPORTED_LOCATION_DEBUG "${COIN_${_libname}_LIBRARY_DEBUG}")
      set_property(TARGET CoinOR::${_libname} PROPERTY INCLUDE_DIRECTORIES "${INCLUDE_DIRS}")
      set_property(TARGET CoinOR::${_libname} PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${INCLUDE_DIRS}")
      # TODO there are probably dependencies across the single libs but we don't care for now. We always included all.
      set_property(TARGET CoinOR::CoinOR APPEND PROPERTY INTERFACE_LINK_LIBRARIES CoinOR::${_libname})
    endif()
  endif()
endmacro()


if(NOT TARGET CoinOR::CoinOR)
  add_library(CoinOR::CoinOR INTERFACE IMPORTED)
  if (VCPKG_TOOLCHAIN)
    # Currently coin-or from vcpkg requires BLAS and LAPACK
    # TODO: Find a better way to do this. Ideal would be if Coin exports a CMake config
    #  Maybe we can parse a header file? Or try_compile?
    #  The current approach fails if VCPKG toolchain is used but CMake somehow finds
    #  an external coin-or. Should be rare to impossible.
    find_package(BLAS)
    find_package(LAPACK)
    target_link_libraries(CoinOR::CoinOR INTERFACE BLAS::BLAS LAPACK::LAPACK)
  endif()
endif()

_coin_find_lib("CBC" "Cbc" "libCbc;Cbc" "libCbcd;Cbc")
_coin_find_lib("CGL" "Cgl" "libCgl;Cgl" "libCgld;Cgl")
_coin_find_lib("CLP" "Clp" "libClp;Clp" "libClpd;Clp")
_coin_find_lib("COINUTILS" "CoinUtils" "libCoinUtils;CoinUtils" "libCoinUtilsd;CoinUtils")
_coin_find_lib("OSI" "Osi" "libOsi;Osi" "libOsid;Osi")
_coin_find_lib("OSI_CLP" "Clp" "libOsiClp;OsiClp" "libOsiClpd;OsiClp")

# TODO allow for COMPONENTS and version parsing/checking
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
  COIN_CBC_LIBRARY
  COIN_CGL_LIBRARY
  COIN_CLP_LIBRARY
  COIN_COINUTILS_LIBRARY
  COIN_OSI_LIBRARY
  COIN_OSI_CLP_LIBRARY
)

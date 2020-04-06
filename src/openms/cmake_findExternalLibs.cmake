# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
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
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# This cmake file handles finding external libs for OpenMS
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# set which library extensions are preferred (we prefer shared libraries)
if(NOT MSVC)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
endif()
if (APPLE)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib;.a")
endif()


#------------------------------------------------------------------------------
# find libs (for linking)
# On Windows:
#   * on windows we need the *.lib versions (dlls alone won't do for linking)
#   * never mix Release/Debug versions of libraries. Leads to strange segfaults,
#     stack corruption etc, due to different runtime libs ...
# compiler-wise: use the same compiler for contrib and OpenMS!
find_package(XercesC REQUIRED)

#------------------------------------------------------------------------------
# BOOST
find_boost(iostreams date_time math_c99 regex)

if(Boost_FOUND)
  message(STATUS "Found Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}" )
  set(CF_OPENMS_BOOST_VERSION_MAJOR ${Boost_MAJOR_VERSION})
	set(CF_OPENMS_BOOST_VERSION_MINOR ${Boost_MINOR_VERSION})
  set(CF_OPENMS_BOOST_VERSION_SUBMINOR ${Boost_SUBMINOR_VERSION})
	set(CF_OPENMS_BOOST_VERSION ${Boost_VERSION})
else()
  message(FATAL_ERROR "Boost or one of its components not found!")
endif()

#------------------------------------------------------------------------------
# SEQAN
find_package(SEQAN 1.4.0 REQUIRED)

#------------------------------------------------------------------------------
# libsvm
find_package(LIBSVM 2.91 REQUIRED) ## will not overwrite LIBSVM_LIBRARY if defined already
if (LIBSVM_FOUND)
	set(CF_OPENMS_LIBSVM_VERSION_MAJOR ${LIBSVM_MAJOR_VERSION})
	set(CF_OPENMS_LIBSVM_VERSION_MINOR ${LIBSVM_MINOR_VERSION})
	set(CF_OPENMS_LIBSVM_VERSION ${LIBSVM_VERSION})
endif()

#------------------------------------------------------------------------------
# COIN-OR
set(CF_USECOINOR 1)
find_package(COIN REQUIRED)

#------------------------------------------------------------------------------
# GLPK
find_package(GLPK REQUIRED)
if (GLPK_FOUND)
	set(CF_OPENMS_GLPK_VERSION_MAJOR ${GLPK_VERSION_MAJOR})
	set(CF_OPENMS_GLPK_VERSION_MINOR ${GLPK_VERSION_MINOR})
	set(CF_OPENMS_GLPK_VERSION ${GLPK_VERSION_STRING})
endif()

#------------------------------------------------------------------------------
# zlib
find_package(ZLIB REQUIRED)

#------------------------------------------------------------------------------
# bzip2
find_package(BZip2 REQUIRED)

#------------------------------------------------------------------------------
# Find eigen3
find_package(Eigen3 3.3.4 REQUIRED)

#------------------------------------------------------------------------------
# Find geometric tools - wildmagick 5
set(WM5_FIND_REQUIRED_COMPONENTS WM5_WM5CORE WM5_WM5MATHEMATICS )
find_package(WM5 REQUIRED)
if (WM5_FOUND)
  add_definitions(${WM5_DEFINITIONS})
endif()

#------------------------------------------------------------------------------
# Find Crawdad libraries if requested
# cmake args: -DCrawdad_DIR=/path/to/Crawdad/ -DWITH_CRAWDAD=TRUE
## TODO check if necessary
if (WITH_CRAWDAD)
  message(STATUS "Will compile with Crawdad support: ${Crawdad_DIR}" )
  find_package(Crawdad REQUIRED)
  # find archive (static) version and add it to the OpenMS library
  find_library(Crawdad_LIBRARY
   NAMES Crawdad.a Crawdad
    HINTS ${Crawdad_DIR})
endif()

#------------------------------------------------------------------------------
# SQLITE
find_package(SQLITE 3.15.0 REQUIRED)

#------------------------------------------------------------------------------
# HDF5
# For MSVC use static linking to the HDF5 libraries
if(MSVC)
  set(HDF5_USE_STATIC_LIBRARIES ON)
endif()
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)

#------------------------------------------------------------------------------
# Done finding contrib libraries
#------------------------------------------------------------------------------

#except for the contrib libs, prefer shared libraries
if(NOT MSVC AND NOT APPLE)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
endif()

#------------------------------------------------------------------------------
# QT
#------------------------------------------------------------------------------
SET(QT_MIN_VERSION "5.5.0")

# find qt
## TODO Use the component variable during install time
## Why were many more QT modules linked? Removed for now until complaints.
set(OpenMS_QT_COMPONENTS Core Network Sql CACHE INTERNAL "QT components for core lib")
find_package(Qt5 COMPONENTS ${OpenMS_QT_COMPONENTS} REQUIRED)

IF (NOT Qt5Core_FOUND)
  message(STATUS "QT5Core not found!")
  message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -DCMAKE_PREFIX_PATH='<path_to_parent_folder_of_lib_folder_withAllQt5Libs>' <src-dir>")
ELSE()
  message(STATUS "Found Qt ${Qt5Core_VERSION}")
ENDIF()


##TODO check if we can integrate the next lines into the openms_add_library cmake macro
add_definitions(${Qt5Core_DEFINITIONS})
add_definitions(${Qt5Network_DEFINITIONS})
add_definitions(${Qt5Sql_DEFINITIONS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Core_EXECUTABLE_COMPILE_FLAGS} ${Qt5Network_EXECUTABLE_COMPILE_FLAGS} ${Qt5Sql_EXECUTABLE_COMPILE_FLAGS}")

# see https://github.com/ethereum/solidity/issues/4124
if("${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}" VERSION_LESS "1.59")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_VARIANT_USE_RELAXED_GET_BY_DEFAULT")
endif()

#------------------------------------------------------------------------------
# PTHREAD
#------------------------------------------------------------------------------
# at least FFSuperHirn requires linking against pthread
find_package (Threads REQUIRED)

#------------------------------------------------------------------------------

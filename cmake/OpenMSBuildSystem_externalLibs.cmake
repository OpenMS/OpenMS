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
# This cmake file handles finding external libs for OpenMS

set(CONTRIB_CUSTOM_DIR CACHE DOC "DEPRECATED: Please use CMAKE_FIND_ROOT_PATH instead! User defined location of contrib dir. If left empty we assume the contrib to be in OpenMS/contrib!")
set(CONTRIB_DIR ${PROJECT_SOURCE_DIR}/contrib/ CACHE INTERNAL "Final contrib path after looking at CMAKE_FIND_ROOT_PATH. Defaults to OpenMS/contrib")

if("${CMAKE_FIND_ROOT_PATH}" STREQUAL "")
	if(NOT "${CONTRIB_CUSTOM_DIR}" STREQUAL "")
		message("Please do no longer use -DCONTRIB_CUSTOM_DIR. This option is deprecated. Please use -DCMAKE_FIND_ROOT_PATH instead.")
		list(INSERT CONTRIB_DIR 0 ${CONTRIB_CUSTOM_DIR})
  else()
		# Append a few unusual default search directories for convenience
		# if no FIND ROOT PATH was specified.
		list(APPEND CONTRIB_DIR /opt/local /usr/local)
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

message(STATUS "CMake find root path: ${CMAKE_FIND_ROOT_PATH}")

set(CONTRIB_INCLUDE_DIR "" CACHE INTERNAL "contrib include dir")
set(CONTRIB_LIB_DIR "" CACHE INTERNAL "contrib lib dir")
foreach(CONTRIB_PATH ${CONTRIB_DIR})
  list(APPEND CONTRIB_INCLUDE_DIR "${CONTRIB_PATH}/include/")
  list(APPEND CONTRIB_LIB_DIR "${CONTRIB_PATH}/lib/")
endforeach()
message(STATUS "Contrib search directories:  ${CONTRIB_DIR}")
message(STATUS "Contrib library directories: ${CONTRIB_LIB_DIR}")
message(STATUS "Contrib include directories: ${CONTRIB_INCLUDE_DIR}")


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

OPENMS_CHECKLIB(CONTRIB_XERCESC "xerces-c_3;xerces-c_static_3;libxerces-c;xerces-c" "xerces-c_3D;xerces-c_static_3D;libxerces-c;xerces-c" "xerces_c")

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
FIND_PACKAGE(SEQAN 1.4.0)
if(SEQAN_FOUND)
  message(STATUS "Found SEQAN version ${SEQAN_VERSION_MAJOR}.${SEQAN_VERSION_MINOR}.${SEQAN_VERSION_PATCH}" )
else()
  message(FATAL_ERROR "SeqAn could not be found. Please install it from www.seqan.de or download and install the OpenMS contrib package.")
endif()

#------------------------------------------------------------------------------
# libsvm
if (WIN32) ## find manually on Windows, as find_package() does not know about debug lib
	OPENMS_CHECKLIB(LIBSVM_LIBRARY "libsvm;svm" "libsvmd;svmd;libsvm;svm" "libSVM")
endif()
find_package(libSVM 2.91) ## will not overwrite LIBSVM_LIBRARY if defined already
if (LIBSVM_FOUND)
	message(STATUS "Found LibSVM version " ${LIBSVM_VERSION})
	set(CF_OPENMS_LIBSVM_VERSION_MAJOR ${LIBSVM_MAJOR_VERSION})
	set(CF_OPENMS_LIBSVM_VERSION_MINOR ${LIBSVM_MINOR_VERSION})
	set(CF_OPENMS_LIBSVM_VERSION ${LIBSVM_VERSION})
	set(DEP_LIBSVM_LIBRARY ${LIBSVM_LIBRARY} ${LIBSVM_LIBRARIES}) # combine for consistent use later
else()
	message(FATAL_ERROR "LibSVM not found!")
endif()

#------------------------------------------------------------------------------
# COIN-OR
if (${USE_COINOR})
	set(CF_USECOINOR 1)
	OPENMS_CHECKLIB(CONTRIB_CBC1 "libCbc;Cbc" "libCbcd;Cbc" "COIN-OR Cbc")
	OPENMS_CHECKLIB(CONTRIB_CBC2 "libCgl;Cgl" "libCgld;Cgl" "COIN-OR Cgl")
	OPENMS_CHECKLIB(CONTRIB_CBC3 "libClp;Clp" "libClpd;Clp" "COIN-OR Clp")
	OPENMS_CHECKLIB(CONTRIB_CBC4 "libCoinUtils;CoinUtils" "libCoinUtilsd;CoinUtils" "COIN-OR Utils")
	OPENMS_CHECKLIB(CONTRIB_CBC5 "libOsi;Osi" "libOsid;Osi" "COIN-OR Osi")
	OPENMS_CHECKLIB(CONTRIB_CBC6 "libOsiClp;OsiClp" "libOsiClpd;OsiClp" "COIN-OR OsiClp")
	set(CONTRIB_CBC ${CONTRIB_CBC1} ${CONTRIB_CBC2} ${CONTRIB_CBC3} ${CONTRIB_CBC4} ${CONTRIB_CBC5} ${CONTRIB_CBC6} )
else()
	set(CF_USECOINOR 0)
	set(CONTRIB_CBC)
endif()

#------------------------------------------------------------------------------
# GLPK
find_package(GLPK REQUIRED)
if (GLPK_FOUND)
	message(STATUS "Found GLPK version " ${GLPK_VERSION_STRING})
	set(CF_OPENMS_GLPK_VERSION_MAJOR ${GLPK_VERSION_MAJOR})
	set(CF_OPENMS_GLPK_VERSION_MINOR ${GLPK_VERSION_MINOR})
	set(CF_OPENMS_GLPK_VERSION ${GLPK_VERSION_STRING})
else()
	message(FATAL_ERROR "GLPK not found!")
endif()

#------------------------------------------------------------------------------
# zlib
find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
  message(STATUS "Found zlib version ${ZLIB_VERSION_STRING}")
else()
  message(FATAL_ERROR "zlib not found!")
endif()

#------------------------------------------------------------------------------
# bzip2
find_package(BZip2 REQUIRED)
if (BZIP2_FOUND)
  message(STATUS "Found bzip2 version ${BZIP2_VERSION_STRING}")
else()
  message(FATAL_ERROR "bzip2 not found!")
endif()

#------------------------------------------------------------------------------
# Find eigen3
FIND_PACKAGE(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
  message(STATUS "Found eigen3 version ${EIGEN3_VERSION}")
else()
  message(FATAL_ERROR "eigen3 not found!")
endif()

#------------------------------------------------------------------------------
# Find geometric tools - wildmagick 5
SET(WM5_FIND_REQUIRED_COMPONENTS WM5_WM5CORE WM5_WM5MATHEMATICS )
find_package(WM5 REQUIRED)
if (WM5_FOUND)
  add_definitions(${WM5_DEFINITIONS})
  message(STATUS "Found WM5")
else()
  message(FATAL_ERROR "WM5 not found!")
endif()

#------------------------------------------------------------------------------
# Done finding contrib libraries
#------------------------------------------------------------------------------

if(MSVC)
	## needed to locate libs (put this above ADD_LIBRARY() - otherwise it will not work)
	link_directories(${CONTRIB_LIB_DIR})
endif()

#except for the contrib libs, prefer shared libraries
if(NOT MSVC AND NOT APPLE)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
endif()

#------------------------------------------------------------------------------
# QT
#------------------------------------------------------------------------------
SET(QT_MIN_VERSION "4.5.0")

# find qt
find_package(Qt4 REQUIRED QtCore QtSql QtNetwork QtGui)

IF (NOT QT4_FOUND)
  message(STATUS "QT4 not found!")
	message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
ENDIF()
include(${QT_USE_FILE})
include(UseQt4)

#------------------------------------------------------------------------------

# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
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
# find libs (for linking)
# On Windows:
#   * on windows we need the *.lib versions (dlls alone won't do for linking)
#   * never mix Release/Debug versions of libraries. Leads to strange segfaults,
#     stack corruption etc, due to different runtime libs ...
# compiler-wise: use the same compiler for contrib and OpenMS!
find_package(XercesC REQUIRED)

#------------------------------------------------------------------------------
# BOOST
set(OpenMS_BOOST_COMPONENTS date_time regex CACHE INTERNAL "Boost components for core lib")
find_boost(iostreams ${OpenMS_BOOST_COMPONENTS})

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
# libsvm
# creates LibSVM target
find_package(LIBSVM 2.91 REQUIRED) ## will not overwrite LIBSVM_LIBRARY if defined already
if (LIBSVM_FOUND)
	set(CF_OPENMS_LIBSVM_VERSION_MAJOR ${LIBSVM_MAJOR_VERSION})
	set(CF_OPENMS_LIBSVM_VERSION_MINOR ${LIBSVM_MINOR_VERSION})
	set(CF_OPENMS_LIBSVM_VERSION ${LIBSVM_VERSION})
endif()

#------------------------------------------------------------------------------
# COIN-OR
# Our find module creates an imported CoinOR::CoinOR target
set(CF_USECOINOR 1)
find_package(COIN REQUIRED)


#------------------------------------------------------------------------------
# GLPK
# creates GLPK target
find_package(GLPK REQUIRED)
if (GLPK_FOUND)
	set(CF_OPENMS_GLPK_VERSION_MAJOR ${GLPK_VERSION_MAJOR})
	set(CF_OPENMS_GLPK_VERSION_MINOR ${GLPK_VERSION_MINOR})
	set(CF_OPENMS_GLPK_VERSION ${GLPK_VERSION_STRING})
endif()

#------------------------------------------------------------------------------
# zlib
# creates ZLIB::ZLIB taret
find_package(ZLIB REQUIRED)

#------------------------------------------------------------------------------
# bzip2
find_package(BZip2 REQUIRED)

#------------------------------------------------------------------------------
# Find eigen3
# creates Eigen3::Eigen3 package
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
# creates SQLite::SQLite3 target
# In our contrib we make a subdir in the includes -> Add PATH_SUFFIXES
# Look for the necessary header
find_path(SQLite3_INCLUDE_DIR NAMES sqlite3.h PATH_SUFFIXES "sqlite")
find_package(SQLite3 3.15.0 REQUIRED)

#------------------------------------------------------------------------------
# HDF5
# For MSVC use static linking to the HDF5 libraries
if(MSVC)
  set(HDF5_USE_STATIC_LIBRARIES ON)
endif()
find_package(HDF5 MODULE REQUIRED COMPONENTS C CXX)

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
SET(QT_MIN_VERSION "5.6.0")

# find qt
set(OpenMS_QT_COMPONENTS Core Network Sql CACHE INTERNAL "QT components for core lib")
find_package(Qt5 ${QT_MIN_VERSION} COMPONENTS ${OpenMS_QT_COMPONENTS} REQUIRED)

IF (NOT Qt5Core_FOUND)
  message(STATUS "QT5Core not found!")
  message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -DCMAKE_PREFIX_PATH='<path_to_parent_folder_of_lib_folder_withAllQt5Libs>' <src-dir>")
ELSE()
  message(STATUS "Found Qt ${Qt5Core_VERSION}")
ENDIF()

# see https://github.com/ethereum/solidity/issues/4124
if("${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}" VERSION_LESS "1.59")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_VARIANT_USE_RELAXED_GET_BY_DEFAULT")
endif()

#------------------------------------------------------------------------------
# PTHREAD
#------------------------------------------------------------------------------
find_package (Threads REQUIRED)


if (WITH_GUI)
  # --------------------------------------------------------------------------
  # Find additional Qt libs
  #---------------------------------------------------------------------------
  set (TEMP_OpenMS_GUI_QT_COMPONENTS Gui Widgets Svg)

  # On macOS the platform plugin of QT requires PrintSupport. We link
  # so it's packaged via the bundling/dependency tools/scripts
  if (APPLE)
    set (TEMP_OpenMS_GUI_QT_COMPONENTS ${TEMP_OpenMS_GUI_QT_COMPONENTS} PrintSupport)
  endif()

  set(OpenMS_GUI_QT_COMPONENTS ${TEMP_OpenMS_GUI_QT_COMPONENTS} CACHE INTERNAL "QT components for GUI lib")

  set(OpenMS_GUI_QT_COMPONENTS_OPT WebEngineWidgets)

  find_package(Qt5 REQUIRED COMPONENTS ${OpenMS_GUI_QT_COMPONENTS})

  IF (NOT Qt5Widgets_FOUND OR NOT Qt5Gui_FOUND OR NOT Qt5Svg_FOUND)
    message(STATUS "Qt5Widgets not found!")
    message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -DCMAKE_PREFIX_PATH='<path_to_parent_folder_of_lib_folder_withAllQt5Libs>' <src-dir>")
  ENDIF()

  find_package(Qt5 QUIET COMPONENTS ${OpenMS_GUI_QT_COMPONENTS_OPT})

  # TODO only works if WebEngineWidgets is the only optional component
  set(OpenMS_GUI_QT_FOUND_COMPONENTS_OPT)
  if(Qt5WebEngineWidgets_FOUND)
    list(APPEND OpenMS_GUI_QT_FOUND_COMPONENTS_OPT "WebEngineWidgets")
  else()
    message(WARNING "Qt5WebEngineWidgets not found, disabling JS Views in TOPPView!")
  endif()
    

  set(OpenMS_GUI_DEP_LIBRARIES "OpenMS")

  foreach(COMP IN LISTS OpenMS_GUI_QT_COMPONENTS)
    list(APPEND OpenMS_GUI_DEP_LIBRARIES "Qt5::${COMP}")
  endforeach()

  foreach(COMP IN LISTS OpenMS_GUI_QT_FOUND_COMPONENTS_OPT)
    list(APPEND OpenMS_GUI_DEP_LIBRARIES "Qt5::${COMP}")
  endforeach()

endif()
#------------------------------------------------------------------------------

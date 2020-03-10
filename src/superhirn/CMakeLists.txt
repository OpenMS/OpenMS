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
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

project("SuperHirn")
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

#------------------------------------------------------------------------------
# Find Boost lib for shared pointers
#------------------------------------------------------------------------------
find_boost()

if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost was not found!")
endif()

# --------------------------------------------------------------------------
# the collection of files for SuperHirn
set(SuperHirn_sources   CACHE INTERNAL "This variable should hold all SuperHirn sources at the end of the config step" )
set(SuperHirn_sources_h CACHE INTERNAL "This variable should hold all SuperHirn headers at the end of the config step" )

# --------------------------------------------------------------------------
# include the file collections
include(source/TRANSFORMATIONS/FEATUREFINDER/sources.cmake)
include(source/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/sources.cmake)
include(include/OpenMS/TRANSFORMATIONS/FEATUREFINDER/sources.cmake)
include(include/OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/sources.cmake)

#------------------------------------------------------------------------------
# QT dependencies
find_package(Qt5 COMPONENTS Core Network REQUIRED)
if (NOT Qt5Network_FOUND)
  message(STATUS "QtNetwork module not found!")
  message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -D QT_QMAKE_EXECUTABLE='<path_to_qmake(.exe)' <src-dir>")
endif()

#------------------------------------------------------------------------------
# add the library
openms_add_library(TARGET_NAME SuperHirn
                   SOURCE_FILES ${SuperHirn_sources}
                   HEADER_FILES ${SuperHirn_sources_h}
                   INTERNAL_INCLUDES ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include
                   EXTERNAL_INCLUDES ${Boost_INCLUDE_DIRS} ${OpenMS_INCLUDE_DIRECTORIES}
                   LINK_LIBRARIES OpenMS
                   DLL_EXPORT_PATH "OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/")

# --------------------------------------------------------------------------
# link libraries for SuperHirn tools
set(SuperHirn_LIBRARIES SuperHirn OpenMS CACHE INTERNAL "SuperHirn libraries" FORCE)

#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
openms_doc_path("${PROJECT_SOURCE_DIR}/include")

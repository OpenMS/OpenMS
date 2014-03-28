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
include(source/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/sources.cmake)
include(source/TRANSFORMATIONS/FEATUREFINDER/sources.cmake)
include(include/OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/sources.cmake)
include(include/OpenMS/TRANSFORMATIONS/FEATUREFINDER/sources.cmake)

# --------------------------------------------------------------------------
# expose SuperHirn public includes (including transitive dependencies)
set(SuperHirn_INTERNAL_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include
                                           ${PROJECT_BINARY_DIR}/include)

set(SuperHirn_EXTERNAL_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS}
                                           ${OpenMS_INCLUDE_DIRECTORIES})

set(SuperHirn_INCLUDE_DIRECTORIES ${SuperHirn_INTERNAL_INCLUDE_DIRECTORIES}
                                  ${SuperHirn_EXTERNAL_INCLUDE_DIRECTORIES}
    CACHE INTERNAL "SuperHirn include directories" FORCE)

# --------------------------------------------------------------------------
# add OpenMS include directories
include_directories(${SuperHirn_INTERNAL_INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SuperHirn_EXTERNAL_INCLUDE_DIRECTORIES})

# --------------------------------------------------------------------------
# add library
add_library(SuperHirn ${SuperHirn_sources})

# --------------------------------------------------------------------------
# autogenerate config.h
set(SuperHirn_CONFIG_H "include/OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnConfig.h")

include(GenerateExportHeader)
generate_export_header(SuperHirn
                      EXPORT_MACRO_NAME SUPERHIRN_DLL_EXPORT
                      EXPORT_FILE_NAME ${SuperHirn_CONFIG_H})

source_group("Header Files\\OpenMS\\" FILES ${SuperHirn_CONFIG_H})
set(SuperHirn_sources_h ${SuperHirn_sources_h} ${SuperHirn_CONFIG_H} CACHE INTERNAL "This variable should hold all SuperHirn headers at the end of the config step")

# we also want to install SuperHirn
install_library(SuperHirn)
install_headers("${SuperHirn_sources_h};${SuperHirn_CONFIG_H}" SuperHirn)

#------------------------------------------------------------------------------
# copy dll to test/doc bin folder on MSVC systems
copy_dll_to_extern_bin(SuperHirn)

# --------------------------------------------------------------------------
# link SuperHirn against OpenMS
target_link_libraries(SuperHirn OpenMS)

# --------------------------------------------------------------------------
# link libraries for SuperHirn tools
set(SuperHirn_LIBRARIES SuperHirn OpenMS CACHE INTERNAL "SuperHirn libraries" FORCE)

#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
openms_doc_path("${PROJECT_SOURCE_DIR}/include")

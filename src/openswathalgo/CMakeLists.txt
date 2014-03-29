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
# $Maintainer: Stephan Aiche, Chris Bielow, Hannes Röst $
# $Authors: Hannes Röst, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenSWATHAlgo")
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

#------------------------------------------------------------------------------
# Find Boost lib
#------------------------------------------------------------------------------
find_boost()

if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost was not found!")
endif()

#------------------------------------------------------------------------------
# Export OpenSWATHAlgo include directories
#------------------------------------------------------------------------------

set(OpenSwathAlgo_EXTERNAL_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS})

set(OpenSwathAlgo_INTERNAL_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include
                                               ${PROJECT_BINARY_DIR}/include)


set(OpenSwathAlgo_INCLUDE_DIRECTORIES ${OpenSwathAlgo_INTERNAL_INCLUDE_DIRECTORIES}
                                      ${OpenSwathAlgo_EXTERNAL_INCLUDE_DIRECTORIES}
    CACHE INTERNAL "OpenSWATHAlgo include directories" FORCE)

#------------------------------------------------------------------------------
# Include directories
include_directories(${OpenSwathAlgo_INTERNAL_INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${OpenSwathAlgo_EXTERNAL_INCLUDE_DIRECTORIES})

#------------------------------------------------------------------------------
# get information regarding the openswath files
include(source/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoFiles.cmake)

#------------------------------------------------------------------------------
# The OpenSwathAlgo library
add_library(OpenSwathAlgo ${OpenSwathAlgoFiles})

#------------------------------------------------------------------------------
# create dll export macros for OpenSwathAlgo
set(OpenSwathAlgo_CONFIG_H "include/OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoConfig.h")
include(GenerateExportHeader)
generate_export_header(OpenSwathAlgo
                      EXPORT_MACRO_NAME OPENSWATHALGO_DLLAPI
                      EXPORT_FILE_NAME ${OpenSwathAlgo_CONFIG_H})

source_group("Header Files\\OpenMS\\ANALYSIS\\OPENSWATHALGO\\" FILES ${OpenSwathAlgo_CONFIG_H})
set(OpenSwathAlgoHeaders ${OpenSwathAlgoHeaders} ${OpenSwathAlgo_CONFIG_H} CACHE INTERNAL "This variable should hold all SuperHirn headers at the end of the config step")

# we also want to install OpenSwathAlgo
install_library(OpenSwathAlgo)
install_headers("${OpenSwathAlgoHeaders};${OpenSwathAlgo_CONFIG_H}" OpenSwathAlgo)

#------------------------------------------------------------------------------
# copy dll to test/doc bin folder on MSVC systems
copy_dll_to_extern_bin(OpenSwathAlgo)

#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
openms_doc_path("${PROJECT_SOURCE_DIR}/include")

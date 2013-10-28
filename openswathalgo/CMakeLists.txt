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
set(Boost_USE_STATIC_LIBS ${BOOST_USE_STATIC})
set(Boost_USE_MULTITHREADED  ON)
set(Boost_USE_STATIC_RUNTIME OFF)
add_definitions(/DBOOST_ALL_NO_LIB) ## disable auto-linking of boost libs (boost tends to guess wrong lib names)
set(Boost_COMPILER "")

# help boost finding it's packages
set(Boost_ADDITIONAL_VERSIONS "1.47.0" "1.48.0" "1.49.0" "1.50.0" "1.51.0" "1.52.0" "1.53.0" "1.54.0")

# 1st attempt does not explicitly requires boost to enable second check (see below)
find_package(Boost 1.42.0)

if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost was not found!")
endif()

#------------------------------------------------------------------------------
# Export OpenSWATHAlgo include directories
#------------------------------------------------------------------------------
set(OPENSWATHALGO_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/include
                                      ${PROJECT_BINARY_DIR}/include
                                      ${Boost_INCLUDE_DIRS}
    CACHE INTERNAL "OpenSWATHAlgo include directories" FORCE)

#------------------------------------------------------------------------------
# Include directories
include_directories(${OPENSWATHALGO_INCLUDE_DIRECTORIES})

#------------------------------------------------------------------------------
# Configure OpenSwathAlgoConfig.h.in file for dll_export macros etc.
# replace any variables in config.h.in with current values
#------------------------------------------------------------------------------
set (CONFIGURED_OPENSWATHALGOCONFIG_H ${PROJECT_BINARY_DIR}/include/OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoConfig.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoConfig.h.in ${CONFIGURED_OPENSWATHALGOCONFIG_H})

#------------------------------------------------------------------------------
# get information regarding the openswath files
include(source/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoFiles.cmake)

#------------------------------------------------------------------------------
# The OpenSWATHAlgo library
add_library(OpenSwathAlgo ${OpenSwathAlgoFiles})

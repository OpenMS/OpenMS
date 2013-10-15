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

cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)
project("OpenMS_host")

#------------------------------------------------------------------------------
# General CMake definitions & helper
#------------------------------------------------------------------------------
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

## CMake sanity check: sometimes CMAKE_SIZEOF_VOID_P just vanishes when
## updating CMake.
if (NOT CMAKE_SIZEOF_VOID_P)
	message(FATAL_ERROR "'CMAKE_SIZEOF_VOID_P' is undefined. Thus you should"
                      " delete CMakeFiles (the directory) and the"
                      " CMakeCache.txt and rerun CMake again! This is some"
                      " weird CMake bug that seems to appear when updating the"
                      " CMake version.")
endif()

#------------------------------------------------------------------------------
# entries meant to be configured using CMake cache
# - do NOT hardcode them here!
# - edit them within CMakeCache.txt using ccmake or cmake -D VAR=VALUE
#------------------------------------------------------------------------------
# CMAKE_FIND_ROOT_PATH
# CMAKE_BUILD_TYPE
# STL_DEBUG
# DB_TEST
# QT_DB_PLUGIN
# MT_CUDA_BUILD_TYPE
#------------------------------------------------------------------------------
option(MT_ENABLE_CUDA "Enable experimental CUDA support" OFF)
option(MT_ENABLE_TBB "Enable Intel TBB support" OFF)
option(MT_ENABLE_OPENMP "Enable OpenMP support" ON)

#------------------------------------------------------------------------------
# OpenMS Version Information
#------------------------------------------------------------------------------

set(OPENMS_PACKAGE_VERSION_MAJOR "1")
set(OPENMS_PACKAGE_VERSION_MINOR "12")
set(OPENMS_PACKAGE_VERSION_PATCH "0")

#------------------------------------------------------------------------------
# Unsupported systems
#------------------------------------------------------------------------------
if (MINGW OR MSYS)
  message(FATAL_ERROR "MSYS and/or MinGW are not supported! Please use a Visual Studio environment! See Windows build instructions for further information!")
endif()

#------------------------------------------------------------------------------
# MacOS X related library extensions
#------------------------------------------------------------------------------
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# due to some awkward updates of cmake in the 2.8.x series we need to add the
  # following
  find_library(CARBON NAMES Carbon)
  find_package(OpenGL)
  set(APPLE_EXTRA_LIBS ${CARBON} ${OPENGL_LIBRARY})
endif()

#------------------------------------------------------------------------------
# Check architecture 32bit vs. 64bit
#------------------------------------------------------------------------------

if (CMAKE_SIZEOF_VOID_P MATCHES "8")
	set(OPENMS_64BIT_ARCHITECTURE 1 CACHE INTERNAL "Architecture-bits")
	message(STATUS "Architecture: 64 bit")
else()
	set(OPENMS_64BIT_ARCHITECTURE 0 CACHE INTERNAL "Architecture-bits")
	message(STATUS "Architecture: 32 bit")
endif()

#------------------------------------------------------------------------------
# Handle build type
#------------------------------------------------------------------------------
# Set default build type (if not set by user on command line)
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

# Force build type into the cache (needs to be set beforehand)
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)

#------------------------------------------------------------------------------
# CMake Utitlities
#------------------------------------------------------------------------------
# include some OpenMS specific macros
include (${PROJECT_SOURCE_DIR}/cmake/OpenMSBuildSystem_macros.cmake)

#------------------------------------------------------------------------------
# Directories
#------------------------------------------------------------------------------

set(OPENMS_HOST_DIRECTORY "${PROJECT_SOURCE_DIR}")
set(OPENMS_HOST_BINARY_DIRECTORY "${PROJECT_BINARY_DIR}")

# Output directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${OPENMS_HOST_BINARY_DIRECTORY}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${OPENMS_HOST_BINARY_DIRECTORY}/bin")
set(OPENMS_BINARY_DIR "${OPENMS_HOST_BINARY_DIRECTORY}/bin")
set(OPENMS_WIN32_DLL_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# (create an empty list if CMAKE_MODULE_PATH does not exist)
if(NOT CMAKE_MODULE_PATH)
  set(CMAKE_MODULE_PATH)
endif()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

#------------------------------------------------------------------------------
# OpenMS core lib
#------------------------------------------------------------------------------
if (MSVC)
	## use OpenMSd.dll in debug mode
	SET(CMAKE_DEBUG_POSTFIX d)
endif()

#------------------------------------------------------------------------------
# OpenMS core lib
#------------------------------------------------------------------------------
add_subdirectory(openms)

#------------------------------------------------------------------------------
# Documentation
#------------------------------------------------------------------------------
add_subdirectory(doc)

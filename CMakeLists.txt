# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#------------------------------------------------------------------------------
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
# -> used from VTK
# Set a consistent MACOSX_RPATH default across all CMake versions.
# When CMake 2.8.12 is required, change this default to 1.
# When CMake 3.0.0 is required, remove this block (see CMP0042).
if (NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()

#------------------------------------------------------------------------------
# Include the standard Dart testing module
# for Nightly Build log
enable_testing()
include(Dart)

#------------------------------------------------------------------------------
# entries meant to be configured using CMake cache
# - do NOT hardcode them here!
# - edit them within CMakeCache.txt using ccmake or cmake -D VAR=VALUE
#------------------------------------------------------------------------------
# CMAKE_PREFIX_PATH
# CMAKE_BUILD_TYPE
# STL_DEBUG
# ADDRESS_SANITIZER
# OPENMS_COVERAGE
#------------------------------------------------------------------------------
option(MT_ENABLE_OPENMP "Enable OpenMP support" ON)
option(BOOST_USE_STATIC "Use Boost static libraries." ON)
option(HAS_XSERVER "Indicates if an X server is available. If set to Off it will disable certain tests and the doc target." ON)
option(ENABLE_TUTORIALS "Indicates whether tutorials should be build. Note that this also depends on the availability of (pdf)latex." ON)
option(WITH_GUI "Build GUI parts of OpenMS (TOPPView&Co). This requires QtGui." ON)

#------------------------------------------------------------------------------
# Extend module path with our modules
# (create an empty list if CMAKE_MODULE_PATH does not exist)
if(NOT CMAKE_MODULE_PATH)
	set(CMAKE_MODULE_PATH)
endif()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")


#------------------------------------------------------------------------------
# OpenMS Version Information
#------------------------------------------------------------------------------

set(OPENMS_PACKAGE_VERSION_MAJOR "2")
set(OPENMS_PACKAGE_VERSION_MINOR "1")
set(OPENMS_PACKAGE_VERSION_PATCH "0")

#------------------------------------------------------------------------------
# retrieve detailed informations on the working copy (git)
include(GetGitRevisionDescription)
option(GIT_TRACKING "Embedding of Git version control information into the OpenMS library (and OpenMSInfo tool). Note that if activated, every branch-switch or commit will trigger re-generation of the build system (MakeFiles, VS Solutions etc)." ON)
if (GIT_TRACKING)
  ## will configure an include file
  git_short_info(OPENMS_GIT_SHORT_REFSPEC OPENMS_GIT_SHORT_SHA1 OPENMS_GIT_LC_DATE)
endif()

message(STATUS "Building OpenMS ${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}")
if(NOT GIT_TRACKING)
  set(OPENMS_GIT_SHORT_SHA1 "disabled")
  set(OPENMS_GIT_SHORT_REFSPEC "disabled")
  message(STATUS "  [CMake is not tracking Git commits and branching. To enable use '-D GIT_TRACKING=ON'.]")
elseif(OPENMS_GIT_SHORT_REFSPEC EQUAL "GIT-NOTFOUND" OR ${OPENMS_GIT_SHORT_REFSPEC} EQUAL "HEAD-HASH-NOTFOUND")
  set(OPENMS_GIT_SHORT_SHA1 "exported")
  set(OPENMS_GIT_SHORT_REFSPEC "exported")
  message(STATUS "  [CMake is not tracking Git commits and branching. Git not found.]")
else()
  # everything found, print some status information
  message(STATUS "  - Repository revision ${OPENMS_GIT_SHORT_SHA1}")
  message(STATUS "  - Repository branch ${OPENMS_GIT_SHORT_REFSPEC}")
  message(STATUS "  - Repository last change date ${OPENMS_GIT_LC_DATE}")
  message(STATUS "  [CMake is tracking Git commits and branching. To disable use '-D GIT_TRACKING=OFF'.]")
endif()

# add some informations used when building packages
include(cmake/package_definitions.cmake)

#------------------------------------------------------------------------------
# Unsupported systems
#------------------------------------------------------------------------------
if (MINGW OR MSYS)
  message(FATAL_ERROR "MSYS and/or MinGW are not supported! Please use a Visual Studio environment! See Windows build instructions for further information!")
endif()

#------------------------------------------------------------------------------
# All the multi-threading stuff (OpenMP, CUDA, TBB)
#------------------------------------------------------------------------------
include(cmake/multithreading.cmake)

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
# Project specific compiler flags
#------------------------------------------------------------------------------
# Fill this with compile flags that external projects should use as well
# for OpenMS internal flags (not promoted to external compiler flags) append to CMAKE_CXX_FLAGS
# see OpenMS/cmake/OpenMSConfig.cmake.in to see how its configured and used (i.e. as OPENMS_ADDCXX_FLAGS)
set(CF_OPENMS_ADDCXX_FLAGS)
include(cmake/compiler_flags.cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CF_OPENMS_ADDCXX_FLAGS}")

#------------------------------------------------------------------------------
# Enable STL debug mode (GCC only)
#------------------------------------------------------------------------------
option(STL_DEBUG "[GCC only] Enable STL-DEBUG mode (very slow)." OFF)
if(STL_DEBUG)
  include(cmake/stl_debug.cmake)
endif()

#------------------------------------------------------------------------------
# Enable AddressSanitizer
#------------------------------------------------------------------------------
option(ADDRESS_SANITIZER "[Clang/GCC only] Enable AddressSanitizer mode (quite slow)." OFF)
if(ADDRESS_SANITIZER)
  include(cmake/AddressSanitizer.cmake)
endif()

#------------------------------------------------------------------------------
# Enable Coverage
#------------------------------------------------------------------------------
option(OPENMS_COVERAGE "[Clang/GCC only] Enable Coverage testing with gcov/lcov. Requires the Debug build type,
 is slower and uses a lot of disk space (>40GB). Enables the target OpenMS_coverage to generate
 a coverage report in html format." OFF)
if(OPENMS_COVERAGE)
  include(cmake/Coverage.cmake)
endif()

#------------------------------------------------------------------------------
# Enable usage statistics
#------------------------------------------------------------------------------
option(ENABLE_USAGE_STATISTICS "Enable collection of usage statistics." ON)
if(ENABLE_USAGE_STATISTICS)
  message(STATUS "Collection of usage statistics and update notifications enabled.")
  message(STATUS "  If you don't want this information to be transmitted to our update sever, you can:")
  message(STATUS "    - Switch the build variable ENABLE_USAGE_STATISTICS to OFF to remove the functionality at build time.")
  message(STATUS "    - Set the environment variable OPENMS_DISABLE_USAGE_STATISTICS to disable the functionality at runtime.")
endif()

#------------------------------------------------------------------------------
# we build shared libraries
set(BUILD_SHARED_LIBS true)

#------------------------------------------------------------------------------
# Host directory for referencing from subprojects
set(OPENMS_HOST_DIRECTORY "${PROJECT_SOURCE_DIR}")
set(OPENMS_HOST_BINARY_DIRECTORY "${PROJECT_BINARY_DIR}")

#------------------------------------------------------------------------------
# CMake Utitlities
#------------------------------------------------------------------------------
# include some OpenMS specific macros
include (${PROJECT_SOURCE_DIR}/cmake/build_system_macros.cmake)
# .. and some to ease the installation process
include (${PROJECT_SOURCE_DIR}/cmake/install_macros.cmake)
# .. and some to ease preparing the documentation builds
include (${PROJECT_SOURCE_DIR}/cmake/doc_macros.cmake)
# .. and some to ease the installation proces
include (${PROJECT_SOURCE_DIR}/cmake/add_library_macros.cmake)
# .. and some to ease the export process
include (${PROJECT_SOURCE_DIR}/cmake/export_macros.cmake)

#------------------------------------------------------------------------------
# Output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${OPENMS_HOST_BINARY_DIRECTORY}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${OPENMS_HOST_BINARY_DIRECTORY}/bin")
set(OPENMS_BINARY_DIR "${OPENMS_HOST_BINARY_DIRECTORY}/bin")
set(OPENMS_WIN32_DLL_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

#------------------------------------------------------------------------------
# define installation subdirectories to allow for custom installations
# note that all those directories are below CMAKE_INSTALL_PREFIX
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_DOC_DIR share/doc CACHE PATH "Installation directory for documentation")
set(INSTALL_SHARE_DIR share/OpenMS CACHE PATH "Installation directory for shared data")

if(WIN32)
  set(DEF_INSTALL_CMAKE_DIR CMake)
  set(DEF_INSTALL_LIB_DIR bin)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/cmake/OpenMS)
  set(DEF_INSTALL_LIB_DIR lib)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")
set(INSTALL_LIB_DIR ${DEF_INSTALL_LIB_DIR} CACHE PATH "Installation directory for libraries")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE DOC SHARE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

#------------------------------------------------------------------------------
# lib naming when building with msvc
if (MSVC)
	## use OpenMSd.dll in debug mode
	SET(CMAKE_DEBUG_POSTFIX d)
endif()

#------------------------------------------------------------------------------
# setup path for external libs (contrib or system)
include(${OPENMS_HOST_DIRECTORY}/cmake/setup_lib_find_paths.cmake)

#------------------------------------------------------------------------------
# Initialize package building environment
set(VALID_PACKAGE_TYPES "none" # the classical build type make/make install
                        "rpm"  # rpm package
                        "deb"  # deb package
                        "dmg") # a drag&drop dmg package for macosx

set(PACKAGE_TYPE "none" CACHE STRING "Package type (internal): ${VALID_PACKAGE_TYPES}")
is_valid_package(${PACKAGE_TYPE})

#------------------------------------------------------------------------------
# The actual OpenMS code (headers and source files)
#------------------------------------------------------------------------------
add_subdirectory(src)

#------------------------------------------------------------------------------
# Documentation
#------------------------------------------------------------------------------
add_subdirectory(doc)

#------------------------------------------------------------------------------
# GKN Package building
#------------------------------------------------------------------------------
# We will use this path for the KNIME packages and the OS X installer
set(SEARCH_ENGINES_DIRECTORY "" CACHE PATH "Directory containing the search engine executables that should be shipped with OpenMS. Note: We expect the layout from the SVN.")

set(ENABLE_PREPARE_KNIME_PACKAGE OFF CACHE BOOL "If enabled, targets to prepare KNIME packages will be generated. Main target will be 'prepare_knime_package'.")
if (ENABLE_PREPARE_KNIME_PACKAGE)
	include(${OPENMS_HOST_DIRECTORY}/cmake/knime_package_support.cmake)
endif()

#------------------------------------------------------------------------------
# Handle export of targets and install them
openms_export_targes()
install_export_targets()

#------------------------------------------------------------------------------
# build settings and configs for external code
#------------------------------------------------------------------------------
# install routines for ..
if("${PACKAGE_TYPE}" STREQUAL "dmg") # .. MacOSX
  include(cmake/package_dragndrop_dmg.cmake)
elseif("${PACKAGE_TYPE}" STREQUAL "rpm") # .. rpms
  include(cmake/package_rpm.cmake)
elseif("${PACKAGE_TYPE}" STREQUAL "deb") # .. deps
  include(cmake/package_deb.cmake)
endif()

#------------------------------------------------------------------------------
# Some final messages for the users
#------------------------------------------------------------------------------
include(cmake/messages.cmake)

message(STATUS "")
message(STATUS "-----------------------------------------------------------------")
message(STATUS "")
message(STATUS "You have successfully configured OpenMS.")
message(STATUS "")
if (MSVC)
  message(STATUS "Execute the 'targets' project to see prominent targets!")
  message(STATUS "For faster loading see the individual solution files located in")
  message(STATUS "the src/ and doc/ folder and their subfolders")
  message(STATUS "  e.g., src/openms/OpenMS.sln")
  message(STATUS "for the solution containing the OpenMS library.")
else()
  message(STATUS "For a full list of make targets execute:")
  message(STATUS "'make targets'")
endif()
message(STATUS "")
message(STATUS "-----------------------------------------------------------------")
message(STATUS "")

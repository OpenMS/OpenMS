# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
# $Maintainer: Julianus Pfeuffer $
# $Authors: Stephan Aiche, Julianus Pfeuffer $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project("OpenMS_host")

# required for cmake < 3.13 if ccache is used
set(CMAKE_AUTOMOC_COMPILER_PREDEFINES OFF)

#------------------------------------------------------------------------------
# General CMake definitions & helper
#------------------------------------------------------------------------------
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

# Set C++ version
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # require at least gcc 4.8 for full C++11 support
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
    message(FATAL_ERROR "GCC version must be at least 4.8!")
  endif()
endif()

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
option(ENABLE_DOCS "Indicates whether documentation should be built." ON)
option(ENABLE_TUTORIALS "Indicates whether tutorials should be build. Note that this also depends on the availability of (pdf)latex." ON)
option(WITH_GUI "Build GUI parts of OpenMS (TOPPView&Co). This requires QtGui." ON)
if(MSVC)
  option(MT_ENABLE_NESTED_OPENMP "Enable nested parallelism." OFF)
else()
  option(MT_ENABLE_NESTED_OPENMP "Enable nested parallelism." ON)
endif()

#------------------------------------------------------------------------------
# Extend module path with our modules
# (create an empty list if CMAKE_MODULE_PATH does not exist)
if(NOT CMAKE_MODULE_PATH)
	set(CMAKE_MODULE_PATH)
endif()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
## CMake looks for the NSIS script in the module path.
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Windows")

#------------------------------------------------------------------------------
# OpenMS Version Information
#------------------------------------------------------------------------------

set(OPENMS_PACKAGE_VERSION_MAJOR "2")
set(OPENMS_PACKAGE_VERSION_MINOR "4")
set(OPENMS_PACKAGE_VERSION_PATCH "0")
set(OPENMS_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}")

message(STATUS "Building OpenMS ${OPENMS_PACKAGE_VERSION}")
#------------------------------------------------------------------------------
# retrieve detailed informations on the working copy (git)
include(GetGitRevisionDescription)
option(GIT_TRACKING "Embedding of Git version control information into the OpenMS library (and OpenMSInfo tool). Note that if activated, every branch-switch or commit will trigger re-generation of the build system (MakeFiles, VS Solutions etc)." ON)
if (GIT_TRACKING)
  ## will configure an include file
  git_short_info(OPENMS_GIT_SHORT_REFSPEC OPENMS_GIT_SHORT_SHA1 OPENMS_GIT_LC_DATE)
endif()

if(NOT GIT_TRACKING)
  set(OPENMS_GIT_SHORT_SHA1 "disabled")
  set(OPENMS_GIT_SHORT_REFSPEC "disabled")
  message(STATUS "  [CMake is not tracking Git commits and branching ('-D GIT_TRACKING=Off')]")
elseif(OPENMS_GIT_SHORT_REFSPEC EQUAL "GIT-NOTFOUND" OR ${OPENMS_GIT_SHORT_REFSPEC} EQUAL "HEAD-HASH-NOTFOUND")
  set(OPENMS_GIT_SHORT_SHA1 "exported")
  set(OPENMS_GIT_SHORT_REFSPEC "exported")
  message(STATUS "  [CMake cannot track Git commits and branching, since Git was not found. ('-D GIT_TRACKING=On')]")
else()
  # everything found, print some status information
  message(STATUS "  - Repository revision ${OPENMS_GIT_SHORT_SHA1}")
  message(STATUS "  - Repository branch ${OPENMS_GIT_SHORT_REFSPEC}")
  message(STATUS "  - Repository last change date ${OPENMS_GIT_LC_DATE}")
  message(STATUS "  [CMake is tracking Git commits and branching ('-D GIT_TRACKING=On')]")
endif()

# Additional versioning information for packaging
# If it is not built from a release branch (or master) it will add the branch name
# and the last change date to everything that can handle strings in the version (eg. installer filenames).
if(${OPENMS_GIT_SHORT_REFSPEC} MATCHES "disabled" OR ${OPENMS_GIT_SHORT_REFSPEC} MATCHES "exported")
  string(TIMESTAMP DATE "%Y%m%d")
  set(OPENMS_PACKAGE_VERSION_PRERELEASE_IDENTIFIER "${OPENMS_GIT_SHORT_REFSPEC}-${DATE}")
  set(OPENMS_PACKAGE_VERSION_FULLSTRING "${OPENMS_PACKAGE_VERSION}-${OPENMS_PACKAGE_VERSION_PRERELEASE_IDENTIFIER}")
elseif(${OPENMS_GIT_SHORT_REFSPEC} MATCHES "release*" OR ${OPENMS_GIT_SHORT_REFSPEC} MATCHES "master")
  set(OPENMS_PACKAGE_VERSION_PRERELEASE_IDENTIFIER "")
  set(OPENMS_PACKAGE_VERSION_FULLSTRING "${OPENMS_PACKAGE_VERSION}")
else()
  string(REPLACE "/" "-" OPENMS_GIT_SHORT_REFSPEC_REPLACED ${OPENMS_GIT_SHORT_REFSPEC})
  string(REPLACE "_" "-" OPENMS_GIT_SHORT_REFSPEC_REPLACED_AGAIN ${OPENMS_GIT_SHORT_REFSPEC_REPLACED})
  set(OPENMS_GIT_LC_DATE_LIST "${OPENMS_GIT_LC_DATE}")
  separate_arguments(OPENMS_GIT_LC_DATE_LIST)
  list(GET OPENMS_GIT_LC_DATE_LIST 0 OPENMS_GIT_LC_DATE_ONLY)
  set(OPENMS_PACKAGE_VERSION_PRERELEASE_IDENTIFIER "${OPENMS_GIT_SHORT_REFSPEC_REPLACED_AGAIN}-${OPENMS_GIT_LC_DATE_ONLY}")
  set(OPENMS_PACKAGE_VERSION_FULLSTRING "${OPENMS_PACKAGE_VERSION}-${OPENMS_PACKAGE_VERSION_PRERELEASE_IDENTIFIER}")
endif()



#------------------------------------------------------------------------------
# Unsupported systems
#------------------------------------------------------------------------------
if (MINGW OR MSYS)
  message(FATAL_ERROR "MSYS and/or MinGW are not supported! Please use a Visual Studio environment! See Windows build instructions for further information!")
endif()

#------------------------------------------------------------------------------
# All the multi-threading stuff (OpenMP)
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
# Enable SIRIUS_test (SiriusAdapter) 
#------------------------------------------------------------------------------
option(ENABLE_SIRIUS_TEST "Enable test for fingerid (SiriusAdapter)" ON)

#------------------------------------------------------------------------------
# Enable CSI:FingerID_test (SiriusAdapter)
#------------------------------------------------------------------------------
option(ENABLE_FINGERID_TEST "Enable test for fingerid (SiriusAdapter)" OFF)

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
option(ENABLE_UPDATE_CHECK "Enable update check." ON)
if(ENABLE_UPDATE_CHECK)
  message(STATUS "Collection of usage statistics and update notifications enabled.")
  message(STATUS "  If you don't want this information to be transmitted to our update sever, you can:")
  message(STATUS "    - Switch the build variable ENABLE_UPDATE_CHECK to OFF to remove the functionality at build time.")
  message(STATUS "    - Set the environment variable OPENMS_DISABLE_UPDATE_CHECK to disable the functionality at runtime.")
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
# .. and some to ease the installation process
include (${PROJECT_SOURCE_DIR}/cmake/add_library_macros.cmake)
# .. and some to ease the export process
include (${PROJECT_SOURCE_DIR}/cmake/export_macros.cmake)


#------------------------------------------------------------------------------
# Initialize package building environment
set(VALID_PACKAGE_TYPES "none" # the classical build type make/make install
                        "rpm"  # rpm package
                        "deb"  # deb package
                        "dmg"  # a drag&drop dmg package for macosx
                        "nsis") # an NSIS installer for Windows

set(PACKAGE_TYPE "none" CACHE STRING "Package type (internal): ${VALID_PACKAGE_TYPES}")
is_valid_package(${PACKAGE_TYPE})

## Warn here because you need the signing identity before the app bundles are defined in the
## src subfolder
if("${PACKAGE_TYPE}" STREQUAL "dmg")
  if (NOT DEFINED SIGNING_IDENTITY OR "${SIGNING_IDENTITY}" STREQUAL "" OR NOT DEFINED SIGNING_EMAIL OR "${SIGNING_EMAIL}" STREQUAL "")
    message(WARNING "Trying to create dmg image without SIGNING_IDENTITY and SIGNING_EMAIL being set. "
                         "This disables codesigning and might lead to problems with Gatekeeper on target machines.")
  else()
    ## CPack has a special variable for the identity
    set(CPACK_BUNDLE_APPLE_CERT_APP "${SIGNING_IDENTITY}")
  endif()
endif()

#------------------------------------------------------------------------------
# Output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${OPENMS_HOST_BINARY_DIRECTORY}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${OPENMS_HOST_BINARY_DIRECTORY}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${OPENMS_HOST_BINARY_DIRECTORY}/lib")
set(OPENMS_BINARY_DIR "${OPENMS_HOST_BINARY_DIRECTORY}/bin")
set(OPENMS_WIN32_DLL_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

#------------------------------------------------------------------------------
# define installation subdirectories to allow for custom installations
# note that all those directories are below CMAKE_INSTALL_PREFIX
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_SHARE_DIR share/OpenMS CACHE PATH "Installation directory for shared data")
set(INSTALL_DOC_DIR share/doc CACHE PATH "Installation directory for documentation")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

if(WIN32)
  set(DEF_INSTALL_CMAKE_DIR CMake)
  set(DEF_INSTALL_LIB_DIR bin)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/cmake/OpenMS)
  set(DEF_INSTALL_LIB_DIR lib)
endif()

set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")
set(INSTALL_LIB_DIR ${DEF_INSTALL_LIB_DIR} CACHE PATH "Installation directory for libraries")
set(INSTALL_PLUGIN_SUBDIR plugins CACHE PATH "Installation subdirectory under the lib dir for QT plugins")
set(INSTALL_PLUGIN_DIR ${INSTALL_LIB_DIR}/${INSTALL_PLUGIN_SUBDIR})
file(RELATIVE_PATH INSTALL_LIB_PATH_REL_TO_BIN /${INSTALL_BIN_DIR} /${INSTALL_LIB_DIR})
set(QT_PLUGIN_PATH_REL_TO_BIN ${INSTALL_LIB_PATH_REL_TO_BIN}/${INSTALL_PLUGIN_SUBDIR})

# Never use RPATH. We have other means on every platform.
# E.g. not yet compatible with our fix_mac_dependencies script.
SET(CMAKE_MACOSX_RPATH FALSE)
# macOS-specific variable used when installing own libraries/bins. Use relative path to bin here.
SET(CMAKE_INSTALL_NAME_DIR "@executable_path/../${INSTALL_LIB_DIR}")

if(UNIX AND NOT APPLE)
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
if (ENABLE_PREPARE_KNIME_PACKAGE)
  ## Since we rely on copying instead of "installing" for KNIME packaging
  ## we need to do it during building already.
  ## In this case, make sure that QT is in a system path if you want to test it locally.
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
else()
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
endif()

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
endif()

#------------------------------------------------------------------------------
# lib naming when building with msvc & convenience location of build.bat
if (MSVC)
	## use OpenMSd.dll in debug mode
	SET(CMAKE_DEBUG_POSTFIX d)
	## copy build.bat to root of binary dir to enable convenient invokation (instead of typing path to source dir all the time)
	if(NOT ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}"))
		file(COPY "${PROJECT_SOURCE_DIR}/tools/build.bat" DESTINATION "${PROJECT_BINARY_DIR}")
	endif()
endif()

#------------------------------------------------------------------------------
# setup path for external libs (contrib or system)
include(${OPENMS_HOST_DIRECTORY}/cmake/setup_lib_find_paths.cmake)


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
# We will use this path for the KNIME packages and the installers
set(SEARCH_ENGINES_DIRECTORY "" CACHE PATH "The location where thirdparty search engines (X!Tandem, OMSSA, MyriMatch, MSGF+, Fido) are located. This directory should have the same structure as the search engine repository at https://github.com/OpenMS/THIRDPARTY after flattening for your platform.")

set(ENABLE_PREPARE_KNIME_PACKAGE OFF CACHE BOOL "If enabled, targets to prepare KNIME packages will be generated. Main target will be 'prepare_knime_package'.")
if (ENABLE_PREPARE_KNIME_PACKAGE)
	include(${OPENMS_HOST_DIRECTORY}/cmake/knime_package_support.cmake)
endif()

#------------------------------------------------------------------------------
# Handle export of targets and install them
openms_export_targets()
install_export_targets()

#------------------------------------------------------------------------------
# Installation and packaging:
#------------------------------------------------------------------------------


# packaging routines:
if(NOT "${PACKAGE_TYPE}" STREQUAL "none")
  include(cmake/package_general.cmake)
  if("${PACKAGE_TYPE}" STREQUAL "dmg") # .. MacOSX
    include(cmake/package_dragndrop_dmg.cmake)
  elseif("${PACKAGE_TYPE}" STREQUAL "rpm") # .. rpms
    include(cmake/package_rpm.cmake)
  elseif("${PACKAGE_TYPE}" STREQUAL "deb") # .. debs
    include(cmake/package_deb.cmake)
  elseif("${PACKAGE_TYPE}" STREQUAL "nsis") # .. Windows nsis
    # For now we use the old NSIS script AS IS
    # i.e. although CMake runs make install etc. NSIS
    # will use everything from the usual build folder like before.
    # TODO automate more and let CMake generate install commands
    # in the script. Unfortunately we have a lot of specific code in the script.
    include(cmake/package_nsis.cmake)
  endif()
  include(CPack)
  include(cmake/package_components.cmake)
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

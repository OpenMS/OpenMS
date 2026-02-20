# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#------------------------------------------------------------------------------
# naming conventions:
#
# prefix a variable with 'CF_' if it is used to configure a file!
# e.g., CF_LibOpenMSExport

# Version should not have unescaped "_" because it is used in the Latex tutorials
set(CF_OPENMS_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION}" CACHE INTERNAL "OpenMS VERSION" FORCE)
set(CF_OPENMS_PACKAGE_VERSION_FULLSTRING "${OPENMS_PACKAGE_VERSION_FULLSTRING}" CACHE INTERNAL "OpenMS version with additional identifier" FORCE)
set(CF_OPENMS_GIT_BRANCH "${OPENMS_GIT_SHORT_REFSPEC}" CACHE INTERNAL "OpenMS Git branch for pyopenms setup.py" FORCE)

#------------------------------------------------------------------------------
# En/disable assertions
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(CF_OPENMS_ASSERTIONS 1)
else()
    set(CF_OPENMS_ASSERTIONS 0)
endif()
set(CF_OPENMS_ASSERTIONS ${CF_OPENMS_ASSERTIONS} CACHE INTERNAL "Enables debug messages (precondition and postconditions are enabled, a bit slower) - this is NOT changing any compiler flags!" FORCE)


#------------------------------------------------------------------------------
# external third-party libs statically linked to OpenMS
#------------------------------------------------------------------------------
# The following library is in our public API, so create export targets
#install(EXPORT EvergreenTarget DESTINATION ${INSTALL_LIB_DIR}/cmake/OpenMS FILE OpenMSTargets.cmake COMPONENT cmake)
add_subdirectory(extern)

#------------------------------------------------------------------------------
# external libs (contrib or system)
#------------------------------------------------------------------------------
#include(${PROJECT_SOURCE_DIR}/cmake_findExternalLibs.cmake)

#------------------------------------------------------------------------------
# configure config.h
#------------------------------------------------------------------------------
include(${PROJECT_SOURCE_DIR}/configh.cmake)

#------------------------------------------------------------------------------
# Build sub-libraries (order matters: Core first, then IO, then Algo)
# Using include() rather than add_subdirectory() so that CMAKE_CURRENT_SOURCE_DIR
# stays as src/openms/, which the existing sources.cmake files depend on for
# their relative paths (e.g. "source/CONCEPT/Constants.cpp").
#------------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/build_core.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/build_io.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/build_algo.cmake)

#------------------------------------------------------------------------------
# Backward-compatible umbrella INTERFACE target
# Consumers can simply link against "OpenMS" and get all three sub-libraries.
#------------------------------------------------------------------------------
add_library(OpenMS INTERFACE)
# All Core→IO back-dependencies have been broken: Core no longer references any
# IO symbols directly.  The linker can process libraries in a single pass.
target_link_libraries(OpenMS INTERFACE OpenMS_Algo OpenMS_IO OpenMS_Core)
target_include_directories(OpenMS INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
  "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>"
)

# Export the umbrella target for find_package(OpenMS) consumers
openms_register_export_target(OpenMS)
install(TARGETS OpenMS EXPORT OpenMSTargets)

# Provide backward-compatible variables for legacy CMake consumers
set(OpenMS_LIBRARIES OpenMS CACHE INTERNAL "OpenMS libraries" FORCE)
set(OpenMS_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_BINARY_DIR}/include
  CACHE INTERNAL "OpenMS include directories" FORCE)

#------------------------------------------------------------------------------
# since the share basically belongs to OpenMS core we control its installation
# here
# Note: that slash(/) is important here, otherwise the whole directory
#       (not its content) will be copied!
#install_directory(${OPENMS_HOST_DIRECTORY}/share/OpenMS/ ${INSTALL_SHARE_DIR} share)
install(DIRECTORY ${OPENMS_HOST_DIRECTORY}/share/OpenMS/
  DESTINATION ${INSTALL_SHARE_DIR}
  COMPONENT share
  FILE_PERMISSIONS      OWNER_WRITE OWNER_READ
                        GROUP_READ
                        WORLD_READ
  DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_EXECUTE GROUP_READ
                        WORLD_EXECUTE WORLD_READ
  REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
  REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
  REGEX "OpenMS/examples" EXCLUDE
)
install_directory(${OPENMS_HOST_DIRECTORY}/share/OpenMS/examples ${INSTALL_SHARE_DIR} examples)

# Devs also need our custom find modules
install_directory(${OPENMS_HOST_DIRECTORY}/cmake/Modules/ ${INSTALL_CMAKE_DIR}/Modules cmake)
file(COPY ${OPENMS_HOST_DIRECTORY}/cmake/Modules DESTINATION ${OPENMS_HOST_BINARY_DIRECTORY})

#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
openms_doc_path("${PROJECT_SOURCE_DIR}/include")

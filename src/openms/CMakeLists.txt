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
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS")
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

set(CMAKE_AUTOMOC ON) # some OpenMS classes make use of Signal/Slot, e.g. SYSTEM/FileWatcher.cpp
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
# third party libs shipped with OpenMS directly
#------------------------------------------------------------------------------
# The following library is in our public API, so create export targets
#install(EXPORT EvergreenTarget DESTINATION ${INSTALL_LIB_DIR}/cmake/OpenMS FILE OpenMSTargets.cmake COMPONENT cmake)
add_subdirectory(thirdparty)

#------------------------------------------------------------------------------
# external libs (contrib or system)
#------------------------------------------------------------------------------
#include(${PROJECT_SOURCE_DIR}/cmake_findExternalLibs.cmake)

#------------------------------------------------------------------------------
# At this point make a summary of where data and doc will be located:
message(STATUS "Info: CF_OPENMS_DATA_PATH: ${CF_OPENMS_DATA_PATH}")
message(STATUS "Info: CF_OPENMS_DOC_PATH: ${CF_OPENMS_DOC_PATH}")

#------------------------------------------------------------------------------
# configure config.h
#------------------------------------------------------------------------------
include(${PROJECT_SOURCE_DIR}/configh.cmake)

#------------------------------------------------------------------------------
# big include file for headers and cpp files, that fills the OpenMS_sources variable
include (${PROJECT_SOURCE_DIR}/includes.cmake)

#------------------------------------------------------------------------------
# all the dependency libraries are linked into libOpenMS.so
#set(OPENMS_DEP_LIBRARIES  ${COIN_LIBRARIES}
#                          ${LIBSVM_LIBRARIES}
#                          ${XercesC_LIBRARY}
#                          ${Boost_LIBRARIES}
#                          ${WM5_LIBRARIES}
#                          ${BZIP2_LIBRARIES}
#                          ${ZLIB_LIBRARIES}
#                          ${SQLite_LIBRARY}
#                          ${HDF5_CXX_LIBRARIES}
#                          ${GLPK_LIBRARIES}
#                          ${CMAKE_DL_LIBS}
#                          ${Qt5Core_LIBRARIES}
#                          ${Qt5Network_LIBRARIES}
#                          ${Qt5Sql_LIBRARIES})

set(OPENMS_DEP_LIBRARIES Evergreen LibSVM::LibSVM XercesC::XercesC Qt5::Core Qt5::Network Qt5::Sql)

set(OPENMS_DEP_PRIVATE_LIBRARIES WM5::WM5 CoinOR::CoinOR Eigen3::Eigen GLPK::GLPK HDF5::HDF5 Boost::iostreams Boost::date_time Boost::boost Boost::regex
  BZip2::BZip2 ZLIB::ZLIB SQLite::SQLite3)

# xerces requires linking against CoreFoundation&CoreServices on macOS
# TODO check if this is still the case
if(APPLE)
  find_library(CoreFoundation_LIBRARY CoreFoundation )
  find_library(CoreServices_LIBRARY CoreServices )
  set(OPENMS_DEP_LIBRARIES ${OPENMS_DEP_LIBRARIES}
                           ${CoreFoundation_LIBRARY}
                           ${CoreServices_LIBRARY})
endif()

if (OPENMP_FOUND)
  list(APPEND OPENMS_DEP_LIBRARIES OpenMP::OpenMP_CXX)
endif()

if (WITH_CRAWDAD) ## TODO check if still necessary
  list(APPEND OPENMS_DEP_LIBRARIES ${Crawdad_LIBRARY})
endif()


## TODO do the following sorting based on the imported targets and their INCLUDE_DIRECTORY property
#set (OPENMS_CONTRIB_INCLUDE_DIRS "")
#if (OPENMS_CONTRIB_LIBS)
#  set (OPENMS_CONTRIB_INCLUDE_DIRS ${OPENMS_CONTRIB_LIBS}/include)
#endif()

## Sort include directories into contrib ones (prefered in order) and system ones.
#foreach (includedir
#         ${SQLite_INCLUDE_DIR}
#         ${HDF5_INCLUDE_DIRS}
#         ${LIBSVM_INCLUDE_DIRS}
#         ${EIGEN3_INCLUDE_DIR}
#         ${WM5_INCLUDE_DIRS}
#         ${CRAWDAD_INCLUDE_DIRS}
#         ${COIN_INCLUDE_DIRS}
#         ${GLPK_INCLUDE_DIRS}
#         ${ZLIB_INCLUDE_DIRS}
#         ${BZIP2_INCLUDE_DIR}
#         ${XercesC_INCLUDE_DIRS}
#         ${Boost_INCLUDE_DIRS})
#  if (${includedir} MATCHES "${OPENMS_CONTRIB_INCLUDE_DIRS}.*")
#    list(APPEND INCLUDES_FROM_CONTRIB ${includedir})
#  else()
#    list(APPEND INCLUDES_FROM_SYSTEM ${includedir})
#  endif()
#endforeach()

openms_add_library(TARGET_NAME  OpenMS
                   SOURCE_FILES  ${OpenMS_sources}
                   HEADER_FILES  ${OpenMS_sources_h}
                                 ${OpenMS_configured_headers}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include
                                     ${CMAKE_CURRENT_BINARY_DIR}/include
                   PRIVATE_INCLUDES ${EOL_BSPLINE_INCLUDE_DIRECTORY}
                                    ${ISOSPEC_INCLUDE_DIRECTORY}
                                    ${JSON_INCLUDE_DIRECTORY}
                   EXTERNAL_INCLUDES
                   LINK_LIBRARIES OpenSwathAlgo
                                  ${OPENMS_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${OPENMS_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "OpenMS/")

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
install_directory(${OPENMS_HOST_DIRECTORY}/share/OpenMS/examples ${INSTALL_SHARE_DIR}/examples examples)

# Devs also need our custom find modules
install_directory(${OPENMS_HOST_DIRECTORY}/cmake/Modules/ ${INSTALL_LIB_DIR}/cmake/OpenMS/Modules cmake)
file(COPY ${OPENMS_HOST_DIRECTORY}/cmake/Modules DESTINATION ${OPENMS_HOST_BINARY_DIRECTORY})

#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
openms_doc_path("${PROJECT_SOURCE_DIR}/include")

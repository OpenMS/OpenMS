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
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS")
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

#------------------------------------------------------------------------------
# naming conventions:
#
# prefix a variable with 'CF_' if it is used to configure a file!
# e.g., CF_LibOpenMSExport
set(CF_OPENMS_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION_MAJOR}.${OPENMS_PACKAGE_VERSION_MINOR}.${OPENMS_PACKAGE_VERSION_PATCH}" CACHE INTERNAL "OpenMS VERSION" FORCE)

#------------------------------------------------------------------------------
# En/disable assertions
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
	set(CF_OPENMS_ASSERTIONS 1)
else()
	set(CF_OPENMS_ASSERTIONS 0)
endif()
set(CF_OPENMS_ASSERTIONS ${CF_OPENMS_ASSERTIONS} CACHE INTERNAL "Enables debug messages (precondition and postconditions are enabled, a bit slower) - this is NOT changing any compiler flags!" FORCE)


#------------------------------------------------------------------------------
# 3rd party libs shipped with OpenMS directly
#------------------------------------------------------------------------------
add_subdirectory(3rdparty)

#------------------------------------------------------------------------------
# external libs (contrib or system)
#------------------------------------------------------------------------------
include(${PROJECT_SOURCE_DIR}/cmake_findExternalLibs.cmake)

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
# all the dependency libraries are linked into libOpenMS.so, except Qt
set(OPENMS_DEP_LIBRARIES  ${COIN_LIBRARIES}
                          ${GSL_LIBRARIES}
                          ${GSL_CBLAS_LIBRARIES}
                          ${LIBSVM_LIBRARIES}
                          ${XercesC_LIBRARY}
                          ${Boost_LIBRARIES}
                          ${WM5_LIBRARIES}
                          ${BZIP2_LIBRARIES}
                          ${ZLIB_LIBRARIES}
                          ${GLPK_LIBRARIES})

# xerces requires linking against CoreFoundation&CoreServices
if(APPLE)
  find_library(CoreFoundation_LIBRARY CoreFoundation )
  find_library(CoreServices_LIBRARY CoreServices )
  set(OPENMS_DEP_LIBRARIES ${OPENMS_DEP_LIBRARIES}
                           ${CoreFoundation_LIBRARY}
                           ${CoreServices_LIBRARY})
endif()

if (MSVC)
	list(APPEND OPENMS_DEP_LIBRARIES opengl32.lib)
endif()
if (WITH_CRAWDAD)
  list(APPEND OPENMS_DEP_LIBRARIES Crawdad)
endif()

openms_add_library(TARGET_NAME  OpenMS
                   SOURCE_FILES  ${OpenMS_sources}
                   HEADER_FILES  ${OpenMS_sources_h}
                                 ${OpenMS_configured_headers}
                   INTERNAL_INCLUDES ${PROJECT_SOURCE_DIR}/include
                                     ${PROJECT_BINARY_DIR}/include
									 PRIVATE_INCLUDES ${EOL_BSPLINE_INCLUDE_DIRECTORY}
                   EXTERNAL_INCLUDES ${QT_INCLUDES}
                                     ${GSL_INCLUDE_DIRS}
                                     ${Boost_INCLUDE_DIRS}
                                     ${SEQAN_INCLUDE_DIRS}
                                     ${LIBSVM_INCLUDE_DIRS}
                                     ${GLPK_INCLUDE_DIRS}
                                     ${ZLIB_INCLUDE_DIRS}
                                     ${BZIP2_INCLUDE_DIR}
                                     ${EIGEN3_INCLUDE_DIR}
                                     ${WM5_INCLUDE_DIRS}
                                     ${OpenSwathAlgo_INCLUDE_DIRECTORIES}
                                     ${CRAWDAD_INCLUDE_DIRS}
                                     ${XercesC_INCLUDE_DIRS}
                                     ${COIN_INCLUDE_DIRS}
                   LINK_LIBRARIES OpenSwathAlgo
                                  ${APPLE_EXTRA_LIBS}
                                  ${QT_LIBRARIES}
                                  ${OPENMS_DEP_LIBRARIES}
                   DLL_EXPORT_PATH "OpenMS/")

#------------------------------------------------------------------------------
# since the share basically belongs to OpenMS core we control its installation
# here
# Note: that slash(/) is important here, otherwise the whole directory
#       (not its content) will be copied!
install_directory(${OPENMS_HOST_DIRECTORY}/share/OpenMS/ ${INSTALL_SHARE_DIR} share)

#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
openms_doc_path("${PROJECT_SOURCE_DIR}/include")

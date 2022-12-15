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
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

################### Download libtorch

# for libraries for linux/CPU only
include(FetchContent)

if (UNIX)
  FetchContent_Declare(
    Torch
    PREFIX Torch
    URL https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.13.0%2Bcpu.zip
  )
endif()

# create Torch properties like, e.g. ${torch_SOURCE_DIR} etc.
FetchContent_GetProperties(Torch)
if(NOT torch_POPULATED)
  FetchContent_Populate(Torch)
endif()

list(APPEND CMAKE_PREFIX_PATH "${torch_SOURCE_DIR}/share/cmake/Torch")
find_package(Torch REQUIRED)

# add Torch library as dependency to OpenMS_ML library
list(APPEND OpenMS_ML_DEP_LIBRARIES torch)


project("OpenMS_ML")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# --------------------------------------------------------------------------
# entry point for include file for headers and cpp files
include (${PROJECT_SOURCE_DIR}/includes.cmake)

#------------------------------------------------------------------------------
# add the library
openms_add_library(TARGET_NAME OpenMS_ML                            # name of the library
                   SOURCE_FILES ${OpenMSML_sources}                 # library source and header files
                   HEADER_FILES ${OpenMSML_sources_h}
                   INTERNAL_INCLUDES ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include
                   EXTERNAL_INCLUDES
                   LINK_LIBRARIES ${OpenMS_ML_DEP_LIBRARIES}        # link dependencies (torch)
                   DLL_EXPORT_PATH "OpenMS/ML/")

### source group definition (probably for VS)
#source_group("Source Files\\OpenMS\\ML" FILES ${sources})
#source_group("Header Files\\OpenMS\\ML" FILES ${header_dataaccess})

# --------------------------------------------------------------------------
# TODO .additional linker flags 
if (MSVC)
  set (ML_lnk_flags "/FORCE:MULTIPLE /INCREMENTAL:NO /ignore:4006 /ignore:4088")
  set_property(TARGET OpenMS_ML APPEND PROPERTY LINK_FLAGS ${ML_lnk_flags}) ## allow multiple definitions of symbols (e.g. from template instanciations or STL-derived classes)
endif()

#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
#openms_doc_path("${PROJECT_SOURCE_DIR}/include")
#openms_doc_path("${PROJECT_SOURCE_DIR}/source/")

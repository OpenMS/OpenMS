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

# required modules
include(CMakeParseArguments)
include(GenerateExportHeader)

# openms_add_library()
# Create an OpenMS like library
function(openms_add_library)
  set(options GENERATE_EXPORT )
  set(oneValueArgs TARGET_NAME DLL_EXPORT_PATH)
  set(multiValueArgs INTERNAL_INCLUDES EXTERNAL_INCLUDES SOURCE_FILES HEADER_FILES)
  cmake_parse_arguments(openms_add_library "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  #------------------------------------------------------------------------------
  # merge into global exported includes
  set(${openms_add_library_TARGET_NAME}_INCLUDE_DIRECTORIES ${openms_add_library_INTERNAL_INCLUDES}
                                                            ${openms_add_library_EXTERNAL_INCLUDES}
      CACHE INTERNAL "${openms_add_library_TARGET_NAME} include directories" FORCE)

  #------------------------------------------------------------------------------
  # Include directories
  include_directories(${openms_add_library_INTERNAL_INCLUDES})
  include_directories(SYSTEM ${openms_add_library_EXTERNAL_INCLUDES})

  #------------------------------------------------------------------------------
  # Add the library
  add_library(${openms_add_library_TARGET_NAME} ${openms_add_library_SOURCE_FILES})

  #------------------------------------------------------------------------------
  # Generate export header if requested
  if(NOT ${openms_add_library_DLL_EXPORT_PATH} STREQUAL "")
    set(_CONFIG_H "include/${openms_add_library_DLL_EXPORT_PATH}/${openms_add_library_TARGET_NAME}Config.h")
    string(TOUPPER ${openms_add_library_TARGET_NAME} _TARGET_UPPER_CASE)
    include(GenerateExportHeader)
    generate_export_header(${openms_add_library_TARGET_NAME}
                          EXPORT_MACRO_NAME ${_TARGET_UPPER_CASE}_DLLAPI
                          EXPORT_FILE_NAME ${_CONFIG_H})

    string(REGEX REPLACE "/" "\\\\" _fixed_path ${openms_add_library_DLL_EXPORT_PATH})

    # add generated header to visual studio
    source_group("Header Files\\${_fixed_path}" FILES ${_CONFIG_H})
  endif()

  #------------------------------------------------------------------------------
  # we also want to install the library
  install_library(${openms_add_library_TARGET_NAME})
  install_headers("${openms_add_library_HEADER_FILES};${_CONFIG_H}" ${openms_add_library_TARGET_NAME})

  #------------------------------------------------------------------------------
  # copy dll to test/doc bin folder on MSVC systems
  copy_dll_to_extern_bin(${openms_add_library_TARGET_NAME})
endfunction()

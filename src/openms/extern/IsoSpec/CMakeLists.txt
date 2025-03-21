# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright OpenMS Inc. -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-present.
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

project("OpenMS-IsoSpec")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)


set(_ISOSPEC_SOURCES
  cwrapper.cpp
  allocator.cpp  
  dirtyAllocator.cpp 
  isoSpec++.cpp 
  isoMath.cpp 
  marginalTrek++.cpp 
  operators.cpp
  element_tables.cpp
  misc.cpp
)

# export include directory
set(ISOSPEC_INCLUDE_DIRECTORY
  ${PROJECT_SOURCE_DIR}
  CACHE
  INTERNAL "isospec include directories" FORCE
)

set(ISOSPEC_SOURCES )
foreach(_file ${_ISOSPEC_SOURCES})
  # append source files
  set(
    ISOSPEC_SOURCES
    ${ISOSPEC_SOURCES}
    ${PROJECT_SOURCE_DIR}/IsoSpec/${_file}
  )
endforeach()

# TODO Instead of including including the complete sources, make a library,
#  link, install and export
# Commented out until used
#add_library(IsoSpec ${ISOSPEC_SOURCES})
#target_include_directories(IsoSpec SYSTEM INTERFACE ${ISOSPEC_INCLUDE_DIRECTORY})

get_filename_component(ISOSPEC_FOLDER "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(APPEND ${TPLICENSEFILE} "\nIsoSpec:\n")
cat(${ISOSPEC_FOLDER}/LICENSE ${TPLICENSEFILE})

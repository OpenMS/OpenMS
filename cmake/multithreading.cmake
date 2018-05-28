# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#------------------------------------------------------------------------------
# This cmake files bundles all the multithreading related stuff of the OpenMS
# build system.

#------------------------------------------------------------------------------
# TBB
#------------------------------------------------------------------------------
set(MT_TBB_INCLUDE_DIR CACHE PATH "Intel Threading Building Blocks 'include' directory.")
set(MT_TBB_LIBRARY_DIR CACHE PATH "Intel Threading Building Blocks libraries directory.")
message(STATUS "Intel TBB: ${MT_ENABLE_TBB}")
if (MT_ENABLE_TBB)
  find_package(TBB)
  if (NOT TBB_FOUND)
    message(FATAL_ERROR "TBB not found but requested.")
  endif()
endif()

if (TBB_FOUND)
  INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
  add_compile_options(-DOPENMS_HAS_TBB)
endif()

#------------------------------------------------------------------------------
# OpenMP
#------------------------------------------------------------------------------
if (MT_ENABLE_OPENMP)
	find_package(OpenMP)
endif()
message(STATUS "OpenMP: ${MT_ENABLE_OPENMP}")

if (OPENMP_FOUND)
  # do NOT use add_definitions() here, because RC.exe on windows will fail
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	if (NOT MSVC)
		set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endif()

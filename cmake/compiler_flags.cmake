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
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# This cmake file handles all the project specific compiler flags

# allow additional custom compile flags on the cmake command line by using -DMY_CXX_FLAGS="-g -D_GLIBCXX_ASSERTIONS ..."
# useful for e.g. Release with debug symbols on gcc/clang
if (MY_CXX_FLAGS)
  message(STATUS "Adding custom compile flags: '${MY_CXX_FLAGS}'!")
  add_compile_options(${MY_CXX_FLAGS})
endif()

########
########    SIMD flags are now set per-target
########
set(x64_CPU "x86|AMD64") ## CMake returns 'x86-64' on Linux and 'AMD64' on Windows..
message(STATUS "Processor is : ${CMAKE_SYSTEM_PROCESSOR}")

# Compiler identification (now used in openms_add_library for target-specific flags)
if (CMAKE_COMPILER_IS_GNUCXX)
  option(ENABLE_GCC_WERROR "Enable -Werror on gcc compilers" OFF)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_COMPILER_IS_CLANG true CACHE INTERNAL "Is CLang compiler (clang++)")
else()
  set(CMAKE_COMPILER_IS_INTELCXX true CACHE INTERNAL "Is Intel C++ compiler (icpc)")
endif()

## -Wconversion flag for GCC
set(CXX_WARN_CONVERSION OFF CACHE BOOL "Enables warnings for type conversion problems (GCC only)")
message(STATUS "Compiler checks for conversion: ${CXX_WARN_CONVERSION}")
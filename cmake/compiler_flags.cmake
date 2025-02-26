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
if (MY_CXX_FLAGS) ## do not change this name! it's used in configh.cmake
  message(STATUS "Custom compile flags: '${MY_CXX_FLAGS}' will be added to targets via target_compiler_flags.cmake")
  # Note: These flags will be added to targets via the helper functions in target_compiler_flags.cmake
endif()

########
########    deal with SSE/AVX flags
########
set(x64_CPU "x86|AMD64") ## CMake returns 'x86-64' on Linux and 'AMD64' on Windows..
message(STATUS "Processor is : ${CMAKE_SYSTEM_PROCESSOR}")
# if we support more ISA's in the future (MIPS, SPARC), then also update OpenMSOSInfo::getActiveSIMDExtensions
if (MSVC)
  ## enable 'AVX' on x86-64, to achive faster base64 en-/decoding via SIMDe
  ## note: MSVC lacks flags for SSE3/SSE4 (only unofficial ones like /d2archSSE42 are available, but SIMDe does not care about them)
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "${x64_CPU}") 
    ## for SIMDe we need to use explicit compiler flags, which in turn define macros (like '#define __AVX__'), which SIMDe will check for and only then create vectorized code
    ## Disabling AVX will actually make the SIMDe code slower compared to the non-SSE version (for Base64 encoding/decoding at least)
    # Note: SIMD flags are now applied per-target in target_compiler_flags.cmake
  endif()
else()  ## GCC/Clang/AppleClang
  ## enable SSE3 on x86, to achive faster base64 en-/decoding
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "${x64_CPU}") 
    # Note: SIMD flags are now applied per-target in target_compiler_flags.cmake
  endif()
endif()
## do nothing for ARM at the moment, since SIMDe will do the right thing upon detecting ARM: https://github.com/simd-everywhere/simde/blob/master/simde/simde-arch.h#L117
## (and it seems that neon instructions compile without error even if no compile flag is given -- as opposed to x64 intrinsics)

####
####  more flags...
####

if (CMAKE_COMPILER_IS_GNUCXX)

  option(ENABLE_GCC_WERROR "Enable -Werror on gcc compilers" OFF)
  if (ENABLE_GCC_WERROR)
    message(STATUS "Enable -Werror for gcc - note that this may not work on all compilers and system settings!")
    # Note: -Werror flag is now applied per-target in target_compiler_flags.cmake
  endif()


  # Recommended setting for eclipse, see http://www.cmake.org/Wiki/CMake:Eclipse
  if (CMAKE_GENERATOR STREQUAL "Eclipse CDT4 - Unix Makefiles")
    # Note: -fmessage-length=0 flag is now applied per-target in target_compiler_flags.cmake
  endif()
  
elseif (MSVC)
	# do not use add_definitions
	# add definitions also lands in stuff like RC_DEFINITION which tend to fail if you use
	# Eclipse CDT 4 - NMAKE generator
	# use set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} ...") instead

	# Note: MSVC-specific flags are now applied per-target in target_compiler_flags.cmake
	# This includes:
	# - /wd4251 /wd4275 (disable dll-interface warning)
	# - /wd4996 (disable deprecated functions warning)
	# - /wd4661 (disable explicit template instantiation request warning)
	# - /wd4503 (disable decorated name length exceeded warning)
	# - /wd4068 (disable unknown pragma warning)
	# - /bigobj (for large object files)
	# - /MP (use multiple CPU cores)
  
	# Note: MSVC-specific definitions are now applied per-target in target_compiler_flags.cmake
	# This includes:
	# - _SCL_SECURE_NO_WARNINGS
	# - _CRT_SECURE_NO_WARNINGS
	# - _CRT_SECURE_NO_DEPRECATE
	# - OPENMS_XERCESDLL (xerces bug workaround)
	# - NOMINMAX (coinor windows.h include bug workaround)

	## hdf5 linkage for windows (in case we want to build dynamically)
	# add_definitions(-DH5_BUILT_AS_DYNAMIC_LIB)

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang") # using regular Clang or AppleClang

  set(CMAKE_COMPILER_IS_CLANG true CACHE INTERNAL "Is CLang compiler (clang++)")
  # Note: Clang-specific flags are now applied per-target in target_compiler_flags.cmake
else()
  set(CMAKE_COMPILER_IS_INTELCXX true CACHE INTERNAL "Is Intel C++ compiler (icpc)")
endif()

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
  CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
  # Note: fPIC flag is now applied per-target in target_compiler_flags.cmake
endif()

## -Wconversion flag for GCC
set(CXX_WARN_CONVERSION OFF CACHE BOOL "Enables warnings for type conversion problems (GCC only)")
if (CXX_WARN_CONVERSION)
  if (CMAKE_COMPILER_IS_GNUCXX)
    # Note: -Wconversion flag is now applied per-target in target_compiler_flags.cmake
  endif()
endif()
message(STATUS "Compiler checks for conversion: ${CXX_WARN_CONVERSION}")
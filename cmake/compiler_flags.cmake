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

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # workaround for MacOS 10.13 and below which does not support std::visit
  # see https://github.com/OpenMS/OpenMS/issues/5714
  add_definitions(-D_LIBCPP_DISABLE_AVAILABILITY)  
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
    add_compile_options(/arch:AVX)
  endif()
else()  ## GCC/Clang/AppleClang
  ## enable SSE3 on x86, to achive faster base64 en-/decoding
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "${x64_CPU}") 
    add_compile_options(-mssse3)
  endif()
endif()
## do nothing for ARM at the moment, since SIMDe will do the right thing upon detecting ARM: https://github.com/simd-everywhere/simde/blob/master/simde/simde-arch.h#L117
## (and it seems that neon instructions compile without error even if no compile flag is given -- as opposed to x64 intrinsics)

####
####  more flags...
####

if (CMAKE_COMPILER_IS_GNUCXX)

  add_compile_options(-Wall -Wextra
    #-fvisibility=hidden # This is now added as a target property for each library.     
    -Wno-unknown-pragmas
    -Wno-long-long 
    -Wno-unknown-pragmas
    -Wno-unused-function
    -Wno-variadic-macros
    )
  
  option(ENABLE_GCC_WERROR "Enable -Werror on gcc compilers" OFF)
  if (ENABLE_GCC_WERROR)
    add_compile_options(-Werror)
    message(STATUS "Enable -Werror for gcc - note that this may not work on all compilers and system settings!")
  endif()


  # Recommended setting for eclipse, see http://www.cmake.org/Wiki/CMake:Eclipse
  if (CMAKE_GENERATOR STREQUAL "Eclipse CDT4 - Unix Makefiles")
    add_compile_options(-fmessage-length=0)
  endif()
  
elseif (MSVC)
	# do not use add_definitions
	# add definitions also lands in stuff like RC_DEFINITION which tend to fail if you use
	# Eclipse CDT 4 - NMAKE generator
	# use set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} ...") instead

	## disable dll-interface warning
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4251 /wd4275")

	## disable deprecated functions warning (e.g. for POSIX functions)
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4996")

	## disable explicit template instantiation request for partially defined classes
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4661")

	## disable warning: "decorated name length exceeded, name was truncated"
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4503")

  ## disable warning: "unknown pragma" (occurs thousands of times for, e.g. '#pragma clang diagnostic ignored "-Wfloat-equal"')
	set(CF_OPENMS_ADDCXX_FLAGS "${CF_OPENMS_ADDCXX_FLAGS} /wd4068")
  
	## don't warn about unchecked std::copy()
	add_definitions(/D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_DEPRECATE)

	## xerces bug workaround
	add_definitions(/DOPENMS_XERCESDLL)
	
	## coinor windows.h include bug workaround
	add_definitions(/DNOMINMAX)

	## hdf5 linkage for windows (in case we want to build dynamically)
	# add_definitions(-DH5_BUILT_AS_DYNAMIC_LIB)

	## some .obj are huge and won't compile in VS2022 debug otherwise:
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")

	## use multiple CPU cores (if available)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang") # using regular Clang or AppleClang

  set(CMAKE_COMPILER_IS_CLANG true CACHE INTERNAL "Is CLang compiler (clang++)")
  ## Clang since v14.0 will use ffp-contract=on by default, i.e. use fused-multiply-add ops, which change results (by usually giving higher precision).
  ## This should not affect us, since we do not use/allow -mfma flags (or -march=native) on x64, thus Clang cannot use these instructions in the first place for runtime code
  ## It may use them for compile time FMA though, even if the host does not have the FMA instruction set - ** its magic**!
  ## Then again, Apple-clang uses FMA on AppleSilicon by default, for compile and runtime FMA, and we need to switch it off:
  add_compile_options(-ffp-contract=off)
  # add clang specific warning levels
  # we should not use -Weverything routinely https://quuxplusone.github.io/blog/2018/12/06/dont-use-weverything/
  add_compile_options(-Wall -Wextra)
  # .. and disable some of the harmless ones
  add_compile_options(
                  -Wno-sign-conversion
                  # These are warnings of low severity, which are disabled
                  # for now until we are down to a reasonable size of warnings.
                  -Wno-long-long
                  -Wno-padded
                  -Wno-global-constructors
                  -Wno-exit-time-destructors
                  -Wno-weak-vtables
                  -Wno-documentation-unknown-command
                  -Wno-undef
                  -Wno-documentation
                  -Wno-source-uses-openmp
                  -Wno-old-style-cast
                  -Wno-c++98-compat
                  -Wno-c++98-compat-pedantic
                  # These are warnings of moderate severity, which are disabled
                  # for now until we are down to a reasonable size of warnings.
                  -Wno-unknown-warning-option
                  -Wno-double-promotion
                  -Wno-unused-template
                  -Wno-conversion
                  -Wno-float-equal
                  -Wno-switch-enum
                  -Wno-missing-prototypes
                  -Wno-missing-variable-declarations
                  -Wno-deprecated
                  -Wno-deprecated-register
                  -Wno-covered-switch-default
                  -Wno-date-time
                  -Wno-missing-noreturn
                  )
else()
  set(CMAKE_COMPILER_IS_INTELCXX true CACHE INTERNAL "Is Intel C++ compiler (icpc)")
endif()

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
  CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
  if (WITH_FPIC)
    add_compile_options(-fPIC)
  endif()
endif()

## -Wconversion flag for GCC
set(CXX_WARN_CONVERSION OFF CACHE BOOL "Enables warnings for type conversion problems (GCC only)")
if (CXX_WARN_CONVERSION)
  if (CMAKE_COMPILER_IS_GNUCXX)
    add_compile_options(-Wconversion)
  endif()
endif()
message(STATUS "Compiler checks for conversion: ${CXX_WARN_CONVERSION}")
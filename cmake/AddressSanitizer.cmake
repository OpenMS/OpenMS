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
# $Maintainer: Hannes Roest $
# $Authors: Hannes Roest, Timo Sachsenberg $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# This cmake file enables the AddressSanitizer and UndefinedBehaviorSanitizer
# see http://clang.llvm.org/docs/AddressSanitizer.html and https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html
#     http://en.wikipedia.org/wiki/AddressSanitizer

function(add_asan_to_target TARGET_NAME_ARG)
  if(ADDRESS_SANITIZER)
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
      # add compiler flag
      if (MSVC)
        message(WARNING "AddressSanitizer can only be enabled for GCC and Clang.")
      else()
        # add AddressSanitizer for compiler and linker
        target_compile_options("${TARGET_NAME_ARG}" 
          PUBLIC 
            -fsanitize=address,undefined
            -fno-sanitize-recover=all
            -fno-sanitize=vptr
            -fno-omit-frame-pointer)
        target_link_options("${TARGET_NAME_ARG}" 
          PUBLIC 
            -fsanitize=address,undefined
            -fno-sanitize-recover=all
            -fno-sanitize=vptr)
        message(STATUS "AddressSanitizer is on.")
      endif()
    else()
      message(WARNING "AddressSanitizer is supported for OpenMS debug mode only.")
      message(WARNING "Build type is ${CMAKE_BUILD_TYPE}")
    endif()
  endif()
endfunction()

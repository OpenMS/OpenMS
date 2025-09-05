# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2023.
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
# $Maintainer: Tom David Müller $
# $Author: Tom David Müller $
# --------------------------------------------------------------------------

file(APPEND ${TPLICENSEFILE} "\nONNX Runtime:\n")
cat(${CMAKE_CURRENT_SOURCE_DIR}/LICENSE ${TPLICENSEFILE})

# Determine platform directory
set(arm_CPU "aarch64|arm64|ARM64")
if(APPLE)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "${arm_CPU}")
        set(ONNXRUNTIME_PLATFORM_DIR "osx-arm64")
    else()
        set(ONNXRUNTIME_PLATFORM_DIR "osx-x86_64") 
    endif()
elseif(WIN32)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "${arm_CPU}")
        set(ONNXRUNTIME_PLATFORM_DIR "win-arm64")
    else()
        set(ONNXRUNTIME_PLATFORM_DIR "win-x64")
    endif()
else()
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "${arm_CPU}")
        set(ONNXRUNTIME_PLATFORM_DIR "linux-aarch64")
    else()
        set(ONNXRUNTIME_PLATFORM_DIR "linux-x64")
    endif()
endif()
message(STATUS "Using ONNX Runtime platform: ${ONNXRUNTIME_PLATFORM_DIR}")

# Determine library
set(ONNXRUNTIME_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/platforms/${ONNXRUNTIME_PLATFORM_DIR}")
if(APPLE)
    set(ONNXRUNTIME_LIBRARY "${ONNXRUNTIME_LIB_DIR}/libonnxruntime.dylib")
elseif(WIN32)
    set(ONNXRUNTIME_LIBRARY "${ONNXRUNTIME_LIB_DIR}/onnxruntime.dll")
else()
    set(ONNXRUNTIME_LIBRARY "${ONNXRUNTIME_LIB_DIR}/libonnxruntime.so")
endif()
message(STATUS "Found ONNX Runtime library: ${ONNXRUNTIME_LIBRARY}")

# Add library
add_library(onnxruntime::onnxruntime SHARED IMPORTED GLOBAL)

# Set library location
set_target_properties(onnxruntime::onnxruntime PROPERTIES
    IMPORTED_LOCATION "${ONNXRUNTIME_LIBRARY}"
)

# On windows we also need the import library
if(WIN32)
    set_target_properties(onnxruntime::onnxruntime PROPERTIES
        IMPORTED_IMPLIB "${ONNXRUNTIME_LIB_DIR}/onnxruntime.lib"
    )
endif()

# Add system include to prevent warnings from ONNX Runtime headers
target_include_directories(onnxruntime::onnxruntime 
    SYSTEM INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

message(STATUS "Created imported target: onnxruntime::onnxruntime")
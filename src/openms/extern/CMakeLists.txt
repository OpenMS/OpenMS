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
# $Maintainer: Chris Bielow $
# $Authors: Stephan Aiche, Chris Bielow $
# --------------------------------------------------------------------------

project("OpenMS-3rd-party-libs")
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

# For use in the subdirs to cat their license
function(cat IN_FILE OUT_FILE)
    file(READ ${IN_FILE} CONTENTS)
    file(APPEND ${OUT_FILE} "${CONTENTS}")
endfunction()

# Prepare a temporary file to "cat" to:
set(TPLICENSEFILE ${OPENMS_HOST_BINARY_DIRECTORY}/thirdparty_licenses.txt)
file(WRITE ${TPLICENSEFILE} "The OpenMS library links statically to the following libraries with the following licenses:\n")
## TODO install it somewhere? Append to main license to show in installer?
##  At least we can use it in conda-forge after configuration

##
## internal packages (not checked if they exist externally)
##
add_subdirectory(evergreen/src) # evergreen Bayesian inference library (0.9)
add_subdirectory(eol-bspline)
add_subdirectory(IsoSpec)
add_subdirectory(GTE)
add_subdirectory(Quadtree)
add_subdirectory(tool_description_lib)

##
## external packages (with fallback option to local version)
##
option(USE_EXTERNAL_JSON "Use an external nlohmann::json library" OFF)
message(STATUS "Option USE_EXTERNAL_JSON: ${USE_EXTERNAL_JSON}")
if(USE_EXTERNAL_JSON)
    find_package(nlohmann_json REQUIRED GLOBAL)  ## global propagates the target up to the parent project
    message(STATUS "External 'NLohmann::json' at: '${nlohmann_json_DIR}'")
else()
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    add_subdirectory(nlohmann_json)
endif()

option(USE_EXTERNAL_SQLITECPP "Use an external SQLiteCpp library" OFF)
message(STATUS "Option USE_EXTERNAL_SQLITECPP: ${USE_EXTERNAL_SQLITECPP}")
if(USE_EXTERNAL_SQLITECPP)
    find_package(SQLiteCpp REQUIRED GLOBAL)      ## global propagates the target up to the parent project
    message(STATUS "External 'SQLiteCpp' at: '${SQLiteCpp_DIR}'")
else()
    add_subdirectory(SQLiteCpp)
endif()

option(USE_EXTERNAL_SIMDE "Use an external SIMDe library" OFF)
message(STATUS "Option USE_EXTERNAL_SIMDE: ${USE_EXTERNAL_SIMDE}")
if(USE_EXTERNAL_SIMDE)
    find_package(SIMDe REQUIRED GLOBAL)  ## global propagates the target up to the parent project
    message(STATUS "External 'SIMDe' at: '${SIMDe_INCLUDE_DIR}'")
else()
    add_subdirectory(simde)
endif()

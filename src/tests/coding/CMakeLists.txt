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

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_coding_tests")

# --------------------------------------------------------------------------
# OpenMS internal style tests
if (NOT WIN32)
  # Check for unintended copies in C++11 for loops
  add_test(
    NAME potential_for_loop_copies
    COMMAND sh -c "grep -R for --exclude-dir=contrib --exclude-dir=cppcheck --exclude-dir=thirdparty ${OPENMS_HOST_DIRECTORY} \
    | grep auto | grep -v begin | grep -v areaBegin | grep -v MZBegin | grep -v OMS_CODING_TEST_EXCLUDE | grep ':' \
    | grep -v '&' | grep -v '//' | grep 'for[ ]*(' || echo 'Success' "
  )

  set_tests_properties(potential_for_loop_copies PROPERTIES
    PASS_REGULAR_EXPRESSION "Success"
  )
endif()
# --------------------------------------------------------------------------
# cppcheck tests
find_package(cppcheck)
if( CPPCHECK_FOUND )
include( ${PROJECT_SOURCE_DIR}/cppcheck.cmake )
endif(CPPCHECK_FOUND)

# --------------------------------------------------------------------------
# cpplint testing
find_program(PYTHON_EXECUTABLE
         python
         DOC "python executable used to perform coding convention test.")
if("${PYTHON_EXECUTABLE}" STREQUAL "PYTHON_EXECUTABLE-NOTFOUND")
  message(STATUS "Info: The program python could not be found. Coding convention check will not be available! Add the location of python(.exe) to your PATH environment variable.")
else()
  include( ${PROJECT_SOURCE_DIR}/createcpplinttests.cmake )
endif("${PYTHON_EXECUTABLE}" STREQUAL "PYTHON_EXECUTABLE-NOTFOUND")

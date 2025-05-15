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
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg, Stephan Aiche $
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
set(DO_NOT_TEST_THESE_FILES_REGEX "(/(moc|ui)_|/MSNumpress.cpp|CachedMzML.cpp|SpectrumAccessOpenMS.cpp|ChromatogramExtractorAlgorithm.cpp|thirdparty)")
set(IGNORE_FILES_IN_BUILD_DIRECTORY "^${PROJECT_BINARY_DIR}")

# --------------------------------------------------------------------------
# add_cpplint_tests : This macro generates cpplint tests for files in the
#                     given directory.
#
# The function searches for all sources files in the given directory and
# and generates a cpplint tests for each individual file.
macro(add_cpplint_tests _directory)
  # find files in _directory
  file(GLOB_RECURSE _source_files
       RELATIVE ${OPENMS_HOST_DIRECTORY}/src/${_directory}/
       ${OPENMS_HOST_DIRECTORY}/src/${_directory}/*.cpp)

  # add tests
  foreach(_file_to_test ${_source_files})
    string( REGEX MATCH ${DO_NOT_TEST_THESE_FILES_REGEX} _do_not_test ${_file_to_test} )
    string( REGEX MATCH ${IGNORE_FILES_IN_BUILD_DIRECTORY} _is_in_bin_dir ${_file_to_test})
    if(NOT _do_not_test AND NOT _is_in_bin_dir)
      set(_test_name "src/${_directory}/${_file_to_test}_cpplint_test")

      add_test(${_test_name}
        "${PYTHON_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/cpplint.py"
        "--verbose=5"
        "--filter=-readability/namespace,-build/namespaces,-whitespace/empty_loop_body,-build/c++11"
        "${OPENMS_HOST_DIRECTORY}/src/${_directory}/${_file_to_test}")

      set_tests_properties(
        ${_test_name}
        PROPERTIES
        FAIL_REGULAR_EXPRESSION
        "${CPPLINT_FAIL_REGULAR_EXPRESSION}")
    endif()
  endforeach()
endmacro()

# --------------------------------------------------------------------------
# create tests for all files in the individual file groups
add_cpplint_tests("openswathalgo")
add_cpplint_tests("openms")
add_cpplint_tests("openms_gui")
add_cpplint_tests("topp")

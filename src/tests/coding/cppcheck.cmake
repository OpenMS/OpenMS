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

include(CppcheckTargets)

## we use only source files for cppcheck
set(SOURCE_FILE_REGEX "\\.cpp$")
# Exclude certain files 
# - MSNumpress.cpp since it is an external file
# - DIAScoring.cpp since it triggers an cppcheck bug (may be fixed in newer versions of cppcheck)
# - MSFraggerAdapter.cpp since it triggers an cppcheck bug (may be fixed in newer versions of cppcheck)
# - third party sources
set(DO_NOT_TEST_THESE_FILES_REGEX "(/(moc|ui)_|/MSNumpress.cpp|/DIAScoring.cpp|MSFraggerAdapter.cpp|/XFDR.cpp|thirdparty)")

# --------------------------------------------------------------------------
# add_cpp_check_tests : This macro generates cppcheck tests for files in the
#                       given directory.
#
# The function searches for all sources files in the given directory and
# and generates a cppcheck tests for each individual file.
macro(add_cpp_check_tests _directory)
  # find files in _directory
  file(GLOB_RECURSE _source_files
       RELATIVE ${OPENMS_HOST_DIRECTORY}/src/${_directory}/
       ${OPENMS_HOST_DIRECTORY}/src/${_directory}/*.cpp)

  # add tests
  foreach(_file_to_test ${_source_files})
    string( REGEX MATCH ${SOURCE_FILE_REGEX} _is_source_file ${_file_to_test} )
    string( REGEX MATCH ${DO_NOT_TEST_THESE_FILES_REGEX} _do_not_test ${_file_to_test} )

    if(_is_source_file AND NOT _do_not_test)
      set(_test_name "src/${_directory}/${_file_to_test}")
      add_cppcheck_sources(${_test_name}
                           ${OPENMS_HOST_DIRECTORY}/src/${_directory}/${_file_to_test}
                           STYLE
                           PERFORMANCE
                           INLINE_SUPPRESSION
                           FAIL_ON_WARNINGS)
    endif()
  endforeach()
endmacro()

# --------------------------------------------------------------------------
# include for cppcheck
set(CPPCHECK_INCLUDEPATH_ARG ${OPENMS_GUI_INCLUDE_DIRS})

# --------------------------------------------------------------------------
add_cpp_check_tests("openswathalgo")
add_cpp_check_tests("openms")
add_cpp_check_tests("openms_gui")
add_cpp_check_tests("topp")

# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry               
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

set(DO_NOT_TEST_THESE_FILES_REGEX "/(moc|ui)_")
set(IGNORE_FILES_IN_BUILD_DIRECTORY "^${PROJECT_BINARY_DIR}")

MACRO(ADD_CPPLINT_TEST FILE_TO_TEST)
  string( REGEX MATCH ${DO_NOT_TEST_THESE_FILES_REGEX} do_not_test ${FILE_TO_TEST} )
  string( REGEX MATCH ${IGNORE_FILES_IN_BUILD_DIRECTORY} is_in_bin_dir ${FILE_TO_TEST})
  if(NOT do_not_test AND NOT is_in_bin_dir)
    add_test(${FILE_TO_TEST}_cpplint_test
      "${PYTHON_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/cmake/cpplint.py"
      "--verbose=5"
      "${PROJECT_SOURCE_DIR}/${FILE_TO_TEST}")

    set_tests_properties(
      ${FILE_TO_TEST}_cpplint_test
      PROPERTIES
      FAIL_REGULAR_EXPRESSION
      "${CPPLINT_FAIL_REGULAR_EXPRESSION}")
  endif()
ENDMACRO()

## create tests for all files in the individual file groups
foreach(i ${OpenMS_sources})
  add_cpplint_test(${i})
endforeach()

foreach(i ${OpenMSVisual_sources})
  add_cpplint_test(${i})
endforeach()

foreach(i ${TOPP_executables})
  add_cpplint_test(${TOPP_DIR}/${i}.C)
endforeach()

foreach(i ${UTILS_executables})
  add_cpplint_test(${UTILS_DIR}/${i}.C)
endforeach()

foreach(i ${GUI_executables})
  add_cpplint_test(${GUI_DIR}/${i}.C)
endforeach()


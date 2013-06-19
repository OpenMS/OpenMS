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
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

include(CppcheckTargets)

## we use only source files for cppcheck
set(SOURCE_FILE_REGEX "\\.C$")

MACRO(ADD_CPPCHECK_TEST FILE_TO_TEST)
  string( REGEX MATCH ${SOURCE_FILE_REGEX} is_source_file ${FILE_TO_TEST} )
  if(is_source_file)
    add_cppcheck_sources(${FILE_TO_TEST} ${FILE_TO_TEST} STYLE FAIL_ON_WARNINGS)
  endif(is_source_file)
ENDMACRO()

set(CPPCHECK_INCLUDEPATH_ARG ${OPENMS_INCLUDE_DIRS})

# library checks
set(OpenMS_cppcheck_sources)
foreach(i ${OpenMS_sources})
  add_cppcheck_test(${i})
endforeach()

# GUI library checks
foreach(i ${OpenMSVisual_sources})
  add_cppcheck_test(${i})
endforeach()

# TOPP checks
foreach(i ${TOPP_executables})
  add_cppcheck_test(${TOPP_DIR}/${i}.C)
endforeach()

# UTILS checks
foreach(i ${UTILS_executables})
  add_cppcheck_test(${UTILS_DIR}/${i}.C)
endforeach()

foreach(i ${GUI_executables})
  add_cppcheck_test(${GUI_DIR}/${i}.C)
endforeach()

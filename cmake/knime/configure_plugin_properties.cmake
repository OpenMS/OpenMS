# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

# include helper functions
include ( ${SCRIPT_DIR}common.cmake )

set(required_variables "SOURCE_PATH;TARGET_PATH;OPENMS_VERSION")
check_variables(required_variables)

find_package(Git)
if(GIT_FOUND)
  file(TO_CMAKE_PATH "${SOURCE_PATH}" _OpenMS_CMAKE_PATH)
  execute_process(COMMAND ${GIT_EXECUTABLE} log -n 1 --simplify-by-decoration --pretty=%ai
        WORKING_DIRECTORY ${_OpenMS_CMAKE_PATH}
        ERROR_QUIET
        OUTPUT_VARIABLE OpenMS_WC_LAST_CHANGED_DATE
        OUTPUT_STRIP_TRAILING_WHITESPACE)

  string(REGEX REPLACE "^([0-9]+)-([0-9]+)-([0-9]+) ([0-9]+):([0-9]+).*"
    "\\1\\2\\3\\4\\5" KNIME_DATE "${OpenMS_WC_LAST_CHANGED_DATE}")
  set(CF_OPENMS_VERSION "${OPENMS_VERSION}.${KNIME_DATE}")
else()
  set(CF_OPENMS_VERSION "${OPENMS_VERSION}")
endif()

# create plugin.properties file
configure_file(${SOURCE_PATH}/cmake/knime/plugin.properties.in
               ${TARGET_PATH}/plugin.properties)

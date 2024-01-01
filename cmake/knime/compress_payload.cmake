# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-.
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
# $Maintainer: Julianus Pfeuffer $
# $Authors: Stephan Aiche, Julianus Pfeuffer $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# include helper functions 
include ( ${SCRIPT_DIR}common.cmake )

set(required_variables "ARCH;PLATFORM;PAYLOAD_FOLDER")
check_variables(required_variables)

set(zip_file "${PAYLOAD_FOLDER}/binaries_${PLATFORM}_${ARCH}.zip")
file(TO_NATIVE_PATH ${zip_file} native_zip)

file(GLOB payload_content "${PAYLOAD_FOLDER}/*")
execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar "cf" "${native_zip}" --format=zip ${payload_content}
    WORKING_DIRECTORY "${PAYLOAD_FOLDER}"
)
# As always, CMake prematurely created an unusable new feature. Let's wait 5 more years until it is usable: https://gitlab.kitware.com/cmake/cmake/-/issues/21653
#file(ARCHIVE_CREATE OUTPUT ${native_zip} PATHS ${payload_content} FORMAT zip)

file(REMOVE_RECURSE ${payload_content})
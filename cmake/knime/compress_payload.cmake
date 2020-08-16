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

set(required_variables "ARCH;PLATFORM;PAYLOAD_FOLDER")
check_variables(required_variables)

find_package(Java REQUIRED)
if (NOT Java_JAR_EXECUTABLE)
	message(FATAL_ERROR "Jar executable not found. It is needed to bundle the sources for the KNIME package. Make sure you have the Java SDK installed and in PATH.")
endif()

set(zip_file "${PAYLOAD_FOLDER}/binaries_${PLATFORM}_${ARCH}.zip")
file(TO_NATIVE_PATH ${zip_file} native_zip)

## For the first file/folder we need to create the zip with jar -c
## Otherwise add files by updating the zip with jar -u
set(INITIALLY_CREATE_ZIP On)
file(GLOB payload_content "${PAYLOAD_FOLDER}/*")
foreach(file ${payload_content})
	string(REPLACE "${PAYLOAD_FOLDER}/" "" trimmed_file ${file})
	
	## This means basically "just zip the files" (TODO we could think of using a lighter zip program instead)
	## Generates no new Manifests
	set(JAR_ARGS "fM") 
	if(NOT INITIALLY_CREATE_ZIP)
		set(JAR_ARGS "u${JAR_ARGS}")
	else()
		set(JAR_ARGS "c${JAR_ARGS}")
		set(INITIALLY_CREATE_ZIP Off)
	endif()

	message(STATUS "${Java_JAR_EXECUTABLE} ${JAR_ARGS} ${native_zip} ${trimmed_file}")
	
	# add to zip file
	execute_process(
		COMMAND "${Java_JAR_EXECUTABLE}" "${JAR_ARGS}" "${native_zip}" "${trimmed_file}"
		WORKING_DIRECTORY ${PAYLOAD_FOLDER}
		RESULT_VARIABLE _result
	)
	
	# remove from fs
	if(IS_DIRECTORY ${file})
		execute_process(
			COMMAND ${CMAKE_COMMAND} -E "remove_directory" "${file}"
			WORKING_DIRECTORY ${PAYLOAD_FOLDER}
			RESULT_VARIABLE _result
		)
	else()
		execute_process(
			COMMAND ${CMAKE_COMMAND} -E "remove" "${file}"
			WORKING_DIRECTORY ${PAYLOAD_FOLDER}
			RESULT_VARIABLE _result
		)
	endif()
endforeach()

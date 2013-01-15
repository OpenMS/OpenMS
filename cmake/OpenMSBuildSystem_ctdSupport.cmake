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

# path were the CTDs will be stored
set(CTD_PATH ${PROJECT_BINARY_DIR}/ctds/descriptors)

# path were the executables can be found
set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR})

# payload paths
set(PAYLOAD_PATH ${PROJECT_BINARY_DIR}/ctds/payload)
set(PAYLOAD_BIN_PATH ${PAYLOAD_PATH}/bin)
set(PAYLOAD_LIB_PATH ${PAYLOAD_PATH}/lib)
set(PAYLOAD_SHARE_PATH ${PAYLOAD_PATH}/share)

# script directory
set(SCRIPT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/knime/)

set(ARCH "")
if(OPENMS_64BIT_ARCHITECTURE)
  set(ARCH "64")
else()
  set(ARCH "32")
endif()

set(PLATFORM "")
if (APPLE)
  set(PLATFORM "mac")
elseif(WIN32)
  set(PLATFORM "win")
else()
  set(PLATFORM "lnx")
endif()

## 
function(remove_parameter_from_ctd toolname param)
  set(FILE_CONTENT "")
  file(READ ${CTD_PATH}/${toolname}.ctd FILE_CONTENT)
  foreach(LINE ${FILE_CONTENT})
    message(STATUS ${FILE})
  endforeach()
endfunction()


# create the target directory
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/ctds/)

# create plugin.properties file
configure_file(${PROJECT_SOURCE_DIR}/cmake/knime/plugin.properties.in 
               ${PROJECT_BINARY_DIR}/ctds/plugin.properties)

# copy the icons
file(COPY        ${PROJECT_SOURCE_DIR}/cmake/knime/icons
     DESTINATION ${PROJECT_BINARY_DIR}/ctds/
     PATTERN ".svn" EXCLUDE)
     
# list of all tools that can generate CTDs
set(CTD_executables ${TOPP_executables} ${UTILS_executables})

# remove tools that do not produce CTDs
list(REMOVE_ITEM CTD_executables PhosphoScoring OpenMSInfo FuzzyDiff GenericWrapper InspectAdapter MascotAdapter PILISIdentification PILISModelCV PILISModelTrainer PILISSpectraGenerator SvmTheoreticalSpectrumGeneratorTrainer OpenSwathMzMLFileCacher PepNovoAdapter)

if(APPLE)
  list(REMOVE_ITEM CTD_executables MyriMatchAdapter)
endif(APPLE)

# pseudo-ctd target
add_custom_target(
  create_ctds
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${CTD_PATH}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CTD_PATH}
  DEPENDS TOPP UTILS
)

# call the tools
foreach(TOOL ${CTD_executables})
	add_custom_command(
		TARGET  create_ctds POST_BUILD
		COMMAND ${TOPP_BIN_PATH}/${TOOL} -write_ctd ${CTD_PATH}
    COMMENT "Creating ctd for ${TOOL}"
	)
endforeach()

# remove those parts of the CTDs we cannot model in KNIME
add_custom_target(
  final_ctds
	COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=OMSSAAdapter -DPARAM=omssa_executable -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
	COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=XTandemAdapter -DPARAM=xtandem_executable -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
	COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=XTandemAdapter -DPARAM=default_input_file -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
	COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=IDPosteriorErrorProbability -DPARAM=output_name -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  DEPENDS create_ctds
)

# create final target that collects all sub-calls
add_custom_target(
	prepare_knime_package
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/mimetypes.xml ${CTD_PATH}
  DEPENDS final_ctds
)

add_custom_target(
  create_payload_share
  # remove old directory 
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${PAYLOAD_SHARE_PATH}
  # create new one and fill with the appropriate content
  COMMAND ${CMAKE_COMMAND} -D make_directory ${PAYLOAD_SHARE_PATH}
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_DIRECTORY=${PAYLOAD_SHARE_PATH} -P ${SCRIPT_DIRECTORY}copy_share.cmake
)

add_custom_target(
  prepare_knime_payload_binaries
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${PAYLOAD_PATH}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_PATH}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_BIN_PATH}
)

# copy the binaries
foreach(TOOL ${CTD_executables})
	add_custom_command(
		TARGET  prepare_knime_payload_binaries POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy ${TOPP_BIN_PATH}/${TOOL} "${PAYLOAD_BIN_PATH}/"
    COMMENT "Assemble binary for ${TOOL}"
	)
endforeach()

add_custom_target(
  prepare_knime_payload_libs
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${PAYLOAD_LIB_PATH}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_LIB_PATH}
  # we need the binaries to determine what libraries we need
  DEPENDS prepare_knime_payload_binaries
)

# assemble the libraries, this differs drastically between the different platforms
if (APPLE) 
  add_custom_command(
    TARGET prepare_knime_payload_libs POST_BUILD
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_dependencies.rb -l ${PAYLOAD_LIB_PATH} -b ${PAYLOAD_BIN_PATH}
    COMMENT "Fix the libraries and binaries (mac style)"
  )
elseif(WIN32)
  # assemble required libraries for win32
else()
  # assemble required libraries for lnx
endif()

# TODO: handle the search engines

# handle the binaries.ini
add_custom_target(
  prepare_knime_payload_ini
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D ARCH=${ARCH} -D PLATFORM=${PLATFORM} -D TARGET_DIR=${PAYLOAD_PATH} -D TEMPLATE_FOLDER=${SCRIPT_DIRECTORY} -P ${SCRIPT_DIRECTORY}copy_binaries_ini.cmake
)
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
set(KNIME_PLUGIN_DIRECTORY ${PROJECT_BINARY_DIR}/ctds CACHE PATH "Directory containing the generated plugin-sources for the OpenMS KNIME package")
set(CTD_PATH ${KNIME_PLUGIN_DIRECTORY}/descriptors)

# path were the executables can be found
set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR})
if(WIN32)
  set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR}/$(ConfigurationName))
endif()

# payload paths
set(PAYLOAD_PATH ${KNIME_PLUGIN_DIRECTORY}/payload)
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
file(MAKE_DIRECTORY ${KNIME_PLUGIN_DIRECTORY})

add_custom_target(
  configure_plugin_properties
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_PATH=${KNIME_PLUGIN_DIRECTORY} -D OPENMS_VERSION=${CF_OPENMS_PACKAGE_VERSION} -P ${SCRIPT_DIRECTORY}configure_plugin_properties.cmake
)

# copy the icons
file(COPY        ${PROJECT_SOURCE_DIR}/cmake/knime/icons
     DESTINATION ${KNIME_PLUGIN_DIRECTORY}
     PATTERN ".svn" EXCLUDE)
     
# list of all tools that can generate CTDs
set(CTD_executables ${TOPP_executables} ${UTILS_executables})

# remove tools that do not produce CTDs
list(REMOVE_ITEM CTD_executables PhosphoScoring OpenMSInfo FuzzyDiff GenericWrapper InspectAdapter MascotAdapter PILISIdentification PILISModelCV PILISModelTrainer PILISSpectraGenerator SvmTheoreticalSpectrumGeneratorTrainer OpenSwathMzMLFileCacher PepNovoAdapter)

# pseudo-ctd target
add_custom_target(
  create_ctds
  # we first create the directory to make sure that the remove command does not fail
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CTD_PATH}  
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${CTD_PATH}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CTD_PATH}
  DEPENDS TOPP UTILS
)

# call the tools
foreach(TOOL ${CTD_executables})
  add_custom_command(
    TARGET  create_ctds POST_BUILD
    COMMAND ${TOPP_BIN_PATH}/${TOOL} -write_ctd ${CTD_PATH}
  )
endforeach()

# remove those parts of the CTDs we cannot model in KNIME
add_custom_target(
  final_ctds
  # OMSSAAdapter
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=OMSSAAdapter -DPARAM=omssa_executable -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  # XTandemAdapter 
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=XTandemAdapter -DPARAM=xtandem_executable -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=XTandemAdapter -DPARAM=default_input_file -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  # IDPosteriorErrorProbability
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=IDPosteriorErrorProbability -DPARAM=output_name -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  # MyriMatchAdapter
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=MyriMatchAdapter -DPARAM=myrimatch_executable -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  DEPENDS create_ctds
)

# create final target that collects all sub-calls
add_custom_target(
  prepare_knime_descriptors
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/mimetypes.xml ${CTD_PATH}
  DEPENDS final_ctds
)

add_custom_target(
  prepare_knime_payload_binaries
  # 1st create the directory to make sure that the remove_directory does not fail
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_PATH}  
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${PAYLOAD_PATH}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_PATH}
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_BIN_PATH}
  DEPENDS TOPP UTILS
)

add_custom_target(
  create_payload_share
  # 1st create the directory to make sure that the remove_directory does not fail
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_SHARE_PATH}
  # remove old directory 
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${PAYLOAD_SHARE_PATH}
  # create new one and fill with the appropriate content
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_SHARE_PATH}
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_DIRECTORY=${PAYLOAD_SHARE_PATH} -P ${SCRIPT_DIRECTORY}copy_share.cmake
  DEPENDS prepare_knime_payload_binaries
)

# copy the binaries
foreach(TOOL ${CTD_executables})
  set(tool_path ${TOPP_BIN_PATH}/${TOOL})
  if(WIN32)
    set(tool_path "${tool_path}.exe")
  endif()
  add_custom_command(
    TARGET  prepare_knime_payload_binaries POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${tool_path} "${PAYLOAD_BIN_PATH}/"
  )
endforeach()

add_custom_target(
  prepare_knime_payload_libs
  # 1st create the directory to make sure that the remove_directory does not fail  
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_LIB_PATH}  
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
  )
elseif(WIN32)
  # assemble required libraries for win32
  # OpenMS, OpenMS_GUI, OpenSWATHAlgo, Qt, xerces
  get_target_property(WIN32_DLLLOCATION OpenMS LOCATION)
  get_filename_component(WIN32_DLLPATH ${WIN32_DLLLOCATION} PATH)
  
  add_custom_command(
    TARGET prepare_knime_payload_libs POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${WIN32_DLLPATH}/OpenMS.dll ${PAYLOAD_LIB_PATH}
    COMMAND ${CMAKE_COMMAND} -E copy ${WIN32_DLLPATH}/OpenMS_GUI.dll ${PAYLOAD_LIB_PATH}
    COMMAND ${CMAKE_COMMAND} -E copy ${WIN32_DLLPATH}/OpenSwathAlgo.dll ${PAYLOAD_LIB_PATH}  
  )

  function(copy_library lib target_path)
    string(REGEX REPLACE "lib$" "dll" target_dll "${lib}")
    file(TO_NATIVE_PATH "${target_dll}" target_native)
    add_custom_command(
      TARGET prepare_knime_payload_libs POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy "${target_native}" "${target_path}"
    )
  endfunction()
  
  set(QT_PAYLOAD_LIBS "QTCORE;QTGUI;QTNETWORK;QTOPENGL;QTSQL;QTSVG;QTWEBKIT;PHONON")
  foreach(QT_PAYLOAD_LIB ${QT_PAYLOAD_LIBS})
    set(target_lib "${QT_${QT_PAYLOAD_LIB}_LIBRARY_RELEASE}")
    copy_library(${target_lib}  ${PAYLOAD_LIB_PATH})
  endforeach()

  # xerces
  set(target_lib)
  get_filename_component(xerces_path ${CONTRIB_XERCESC_OPT} PATH)
  file(TO_NATIVE_PATH "${xerces_path}/xerces-c_3_0.dll" target_native)
      add_custom_command(
      TARGET prepare_knime_payload_libs POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy "${target_native}" "${PAYLOAD_LIB_PATH}"
    )
else()
  # assemble required libraries for lnx
	set(QT_PAYLOAD_LIBS "QTCORE;QTGUI;QTNETWORK;QTOPENGL;QTSQL;QTSVG;QTWEBKIT;PHONON")
	foreach(QT_PAYLOAD_LIB ${QT_PAYLOAD_LIBS})
		if(NOT "${QT_${QT_PAYLOAD_LIB}_LIBRARY_RELEASE}" STREQUAL "QT_${QT_PAYLOAD_LIB}_LIBRARY_RELEASE-NOTFOUND")
			set(target_lib "${QT_${QT_PAYLOAD_LIB}_LIBRARY_RELEASE}")
			add_custom_command(
				TARGET prepare_knime_payload_libs POST_BUILD
				COMMAND ${CMAKE_COMMAND} -E copy "${target_lib}" "${PAYLOAD_LIB_PATH}"
			)
		endif()
	endforeach()

	add_custom_command(
		TARGET prepare_knime_payload_libs POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/lib/libOpenMS.so ${PAYLOAD_LIB_PATH}
		COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/lib/libOpenMS_GUI.so ${PAYLOAD_LIB_PATH}
		COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/lib/libOpenSwathAlgo.so ${PAYLOAD_LIB_PATH}  
	)
endif()

# handle the binaries.ini
add_custom_target(
  prepare_knime_payload_ini
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D ARCH=${ARCH} -D PLATFORM=${PLATFORM} -D TARGET_DIR=${PAYLOAD_PATH} -D TEMPLATE_FOLDER=${SCRIPT_DIRECTORY} -P ${SCRIPT_DIRECTORY}copy_binaries_ini.cmake
  DEPENDS prepare_knime_payload_binaries
)

# check if we have valid search engines
if(NOT EXISTS ${SEARCH_ENGINES_DIRECTORY})
  message(FATAL_ERROR "Please specify the path to the search engines to build the KNIME packages. Call cmake -D SEARCH_ENGINES_DIRECTORY=<Path-To-Checkedout-SE> .")
elseif(NOT EXISTS ${SEARCH_ENGINES_DIRECTORY}/OMSSA OR NOT EXISTS ${SEARCH_ENGINES_DIRECTORY}/OMSSA)
  message(FATAL_ERROR "The given search engine directory seems to have an invalid layout. Please check use the one from the SVN.")
elseif(NOT APPLE AND NOT EXISTS ${SEARCH_ENGINES_DIRECTORY}/MyriMatch)
  message(FATAL_ERROR "The given search engine directory seems to have an invalid layout. Please check use the one from the SVN.")
endif()

add_custom_target(
  prepare_knime_payload_searchengines
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SE_PATH=${SEARCH_ENGINES_DIRECTORY} -D TARGET_DIRECTORY=${PAYLOAD_BIN_PATH} -P ${SCRIPT_DIRECTORY}copy_searchengines.cmake
  # We need the folder layout from the bin target
  DEPENDS prepare_knime_payload_binaries
)
  
# the complete payload target
add_custom_target(
  prepare_knime_payload
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D ARCH=${ARCH} -D PLATFORM=${PLATFORM} -D PAYLOAD_FOLDER=${PAYLOAD_PATH} -P ${SCRIPT_DIRECTORY}compress_payload.cmake
  DEPENDS prepare_knime_payload_binaries prepare_knime_payload_libs create_payload_share prepare_knime_payload_ini prepare_knime_payload_searchengines
  )

add_custom_target(
  prepare_meta_information
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/LICENSE ${KNIME_PLUGIN_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/COPYRIGHT ${KNIME_PLUGIN_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/DESCRIPTION ${KNIME_PLUGIN_DIRECTORY}/  
)

add_custom_target(
  prepare_knime_package
  DEPENDS prepare_meta_information configure_plugin_properties prepare_knime_descriptors prepare_knime_payload
  )

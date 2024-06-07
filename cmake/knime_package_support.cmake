# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2023.
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
# $Maintainer: Julianus Pfeuffer$
# $Authors: Stephan Aiche, Julianus Pfeuffer$
# --------------------------------------------------------------------------

# path were the CTDs will be stored
set(KNIME_PROJECT_DIRECTORY ${PROJECT_BINARY_DIR}/knime CACHE PATH "Directory containing the generated plugin-sources for the OpenMS KNIME package")

set(OPENMS_FEATURE_DIRECTORY ${KNIME_PROJECT_DIRECTORY}/de.openms.feature)
set(OPENMS_TP_FEATURE_DIRECTORY ${KNIME_PROJECT_DIRECTORY}/de.openms.thirdparty.feature)

set(KNIME_PLUGIN_DIRECTORY ${OPENMS_FEATURE_DIRECTORY}/de.openms.topp)
set(KNIME_LIB_DIRECTORY ${OPENMS_FEATURE_DIRECTORY}/de.openms.lib)
set(KNIME_TP_PLUGIN_DIRECTORY ${OPENMS_TP_FEATURE_DIRECTORY}/de.openms.thirdparty)
set(TOP_LEVEL_DIRS ${OPENMS_FEATURE_DIRECTORY} ${OPENMS_TP_FEATURE_DIRECTORY})

# create the target directory at configure time
file(MAKE_DIRECTORY ${KNIME_PROJECT_DIRECTORY})

# payload paths
# for the topp plugin
set(CTD_PATH ${KNIME_PLUGIN_DIRECTORY}/descriptors)
set(PAYLOAD_PATH ${KNIME_PLUGIN_DIRECTORY}/payload)
set(PAYLOAD_BIN_PATH ${PAYLOAD_PATH}/bin)

# for the lib plugin
set(LIB_PATH ${KNIME_LIB_DIRECTORY}/payload)
set(PAYLOAD_LIB_PATH ${LIB_PATH}/lib)
set(PAYLOAD_SHARE_PATH ${LIB_PATH}/share)

# for the thirdparty plugin
set(CTD_TP_PATH ${KNIME_TP_PLUGIN_DIRECTORY}/descriptors)
set(TP_PAYLOAD_PATH ${KNIME_TP_PLUGIN_DIRECTORY}/payload)
set(TP_PAYLOAD_BIN_PATH ${TP_PAYLOAD_PATH}/bin)

set(SUB_LEVEL_DIRS ${KNIME_PLUGIN_DIRECTORY} ${KNIME_LIB_DIRECTORY} ${KNIME_TP_PLUGIN_DIRECTORY}
    ${CTD_PATH} ${PAYLOAD_PATH} ${PAYLOAD_BIN_PATH}
    ${LIB_PATH} ${PAYLOAD_LIB_PATH} ${PAYLOAD_SHARE_PATH}
    ${CTD_TP_PATH} ${TP_PAYLOAD_PATH} ${TP_PAYLOAD_BIN_PATH})

# path where the executables can be found
if(CMAKE_CONFIGURATION_TYPES)
  set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR}/${CMAKE_CFG_INTDIR})
else()
  set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR})
endif()

# Use the builtin CMake module to gather paths for Windows runtimes to ship
# Do not install. We rely on copying at build time in this script
if(DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP)
  set(OLD_CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP})
endif()

set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP TRUE)

# collect compiler-provided system runtime libraries (e.g., VS runtime libraries)
include(InstallRequiredSystemLibraries)

# Reset var, just in case
if(DEFINED OLD_CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP)
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP ${OLD_CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP})
endif()

# Find Qt6 includes for KNIME packaging
find_package(Qt6 COMPONENTS ${OpenMS_QT_COMPONENTS} REQUIRED)
get_target_property(QT_QMAKE_EXECUTABLE Qt6::qmake IMPORTED_LOCATION)
exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-query QT_INSTALL_LIBS" OUTPUT_VARIABLE QT_INSTALL_LIBS)
exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-query QT_INSTALL_BINS" OUTPUT_VARIABLE QT_INSTALL_BINS)

# script directory
set(SCRIPT_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/knime/)

# variables for the scripts
set(ARCH "")
if(OPENMS_64BIT_ARCHITECTURE)
  set(ARCH "64")
else()
  set(ARCH "32")
endif()

if(APPLE)
	set(MACOS_TARGET_ARCHS ${CMAKE_OSX_ARCHITECTURES})
	if (NOT MACOS_TARGET_ARCHS)
		# Warning: if cmake is a subprocess of a process that is run under Rosetta,
  		# it will yield x86_64 (but probably also build for it. Therefore it should be fine.)
		set(MACOS_TARGET_ARCHS ${CMAKE_HOST_SYSTEM_PROCESSOR})
	endif()
  # Name according to GenericKNIMENodes specification
  if (MACOS_TARGET_ARCHS STREQUAL "x86_64")
    set(ARCH "64")
  elseif (MACOS_TARGET_ARCHS STREQUAL "arm64")
    set(ARCH "arm64")
  elseif ("x86_64" IN_LIST MACOS_TARGET_ARCHS AND "arm64" IN_LIST MACOS_TARGET_ARCHS)
    set(ARCH "universal")
  else ()
    message(ERROR "Couldn't determine MACOS_TARGET_ARCHS.")
  endif()
endif()

set(PLATFORM "")
if (APPLE)
  set(PLATFORM "mac")
elseif(WIN32)
  set(PLATFORM "win")
else()
  set(PLATFORM "lnx")
endif()

# pseudo-ctd target
add_custom_target(
    create_knime_folders
    DEPENDS TOPP
)

foreach (PATH IN LISTS TOP_LEVEL_DIRS)
  # we first create the directory to make sure that the remove command does not fail
  add_custom_command(
      TARGET  create_knime_folders POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory ${PATH}
      COMMAND ${CMAKE_COMMAND} -E remove_directory ${PATH}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${PATH}
  )
endforeach()

foreach (PATH IN LISTS SUB_LEVEL_DIRS)
  # we know the dirs are empty now
  add_custom_command(
      TARGET  create_knime_folders POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory ${PATH}
  )
endforeach()

find_package(Git)
if(GIT_FOUND)
  file(TO_CMAKE_PATH "${CMAKE_CURRENT_LIST_DIR}" DIR)
  execute_process(COMMAND ${GIT_EXECUTABLE} log -n 1 --simplify-by-decoration --pretty=%ai
                  WORKING_DIRECTORY ${DIR}
                  ERROR_QUIET
                  OUTPUT_VARIABLE OpenMS_WC_LAST_CHANGED_DATE
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  
  string(REGEX REPLACE "^([0-9]+)-([0-9]+)-([0-9]+) ([0-9]+):([0-9]+).*"
         "\\1\\2\\3\\4\\5" KNIME_DATE "${OpenMS_WC_LAST_CHANGED_DATE}")
  set(KNIME_OPENMS_VERSION "${OPENMS_PACKAGE_VERSION}.${KNIME_DATE}")
else()
  set(KNIME_OPENMS_VERSION "${OPENMS_PACKAGE_VERSION}")
endif()

add_custom_target(
  configure_plugin_properties
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_PATH=${KNIME_PLUGIN_DIRECTORY} -D VERSION=${KNIME_OPENMS_VERSION} -D TYPE="TOPPPLUGIN" -P ${SCRIPT_DIRECTORY}configure_plugin_properties.cmake
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_PATH=${KNIME_TP_PLUGIN_DIRECTORY} -D VERSION=${KNIME_OPENMS_VERSION} -D TYPE="TPPLUGIN" -P ${SCRIPT_DIRECTORY}configure_plugin_properties.cmake
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_PATH=${KNIME_LIB_DIRECTORY} -D VERSION=${KNIME_OPENMS_VERSION} -D TYPE="LIBPLUGIN" -P ${SCRIPT_DIRECTORY}configure_plugin_properties.cmake
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_PATH=${OPENMS_FEATURE_DIRECTORY} -D VERSION=${KNIME_OPENMS_VERSION} -D TYPE="FEATURE" -P ${SCRIPT_DIRECTORY}configure_plugin_properties.cmake
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_PATH=${OPENMS_TP_FEATURE_DIRECTORY} -D VERSION=${KNIME_OPENMS_VERSION} -D TYPE="TPFEATURE" -P ${SCRIPT_DIRECTORY}configure_plugin_properties.cmake
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_PATH=${KNIME_PROJECT_DIRECTORY} -D VERSION=${KNIME_OPENMS_VERSION} -D TYPE="UPDATESITE" -P ${SCRIPT_DIRECTORY}configure_plugin_properties.cmake
  DEPENDS create_knime_folders
)

# list of all tools that can generate CTDs and do not include GUI libraries
# TODO make a new category for Adapters
set(CTD_executables ${TOPP_TOOLS})

# remove tools that do not produce CTDs or should not be shipped (because of dependencies or specifics that can not be resolved in KNIME)
list(REMOVE_ITEM CTD_executables OpenMSInfo Resampler ExecutePipeline INIUpdater ImageCreator GenericWrapper MascotAdapter OpenSwathMzMLFileCacher)

# TODO do regex with "Adapter". Safe enough?
set(THIRDPARTY_ADAPTERS
    "MaRaClusterAdapter"
    "XTandemAdapter"
    "MSGFPlusAdapter"
    "LuciphorAdapter"
    "CometAdapter"
    "PercolatorAdapter"
    "MSFraggerAdapter"
    "NovorAdapter"
  )

add_custom_target(
    create_ctds
    DEPENDS create_knime_folders
)

# call the tools to write ctds
foreach(TOOL ${CTD_executables})
  if(${TOOL} IN_LIST THIRDPARTY_ADAPTERS)
    add_custom_command(
      TARGET  create_ctds POST_BUILD
      COMMAND ${TOPP_BIN_PATH}/${TOOL} -write_ctd ${CTD_TP_PATH}
    )
  else()
    add_custom_command(
        TARGET  create_ctds POST_BUILD
        COMMAND ${TOPP_BIN_PATH}/${TOOL} -write_ctd ${CTD_PATH}
    )
  endif()
endforeach()

# Note: We expose FileConverter twice.
# Once as FileConverter in the OpenMS (core) plugin without raw file support (see first if-case in foreach above)
# Once in the Thirdparty plugin as RawFileConverter with raw file support (see below).
# We rename the filename to show a different node name in KNIME but leave the tool name inside the CTD unchanged, so it finds the tool binary.
# TODO change description and accepting file types?
add_custom_command(
  TARGET  create_ctds POST_BUILD
  COMMAND ${TOPP_BIN_PATH}/FileConverter -write_ctd ${CTD_TP_PATH}
  COMMAND ${CMAKE_COMMAND} -E rename ${CTD_TP_PATH}/FileConverter.ctd ${CTD_TP_PATH}/RawFileConverter.ctd
  COMMAND ${CMAKE_COMMAND} -DSCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=RawFileConverter -DCTD_FILE=${CTD_TP_PATH}/RawFileConverter.ctd -P ${SCRIPT_DIRECTORY}change_exec_name_in_ctd.cmake
)

# remove those parts of the CTDs we cannot or do not want to model in KNIME
# e.g. paths to executables that we ship and whose directories are in path environment
add_custom_target(
  final_ctds
  # MaRaClusterAdapter
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=MaRaClusterAdapter -DPARAM=maracluster_executable -D CTD_PATH=${CTD_TP_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  # XTandemAdapter
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=XTandemAdapter -DPARAM=xtandem_executable -D CTD_PATH=${CTD_TP_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  # MSGFPlusAdapter
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=MSGFPlusAdapter -DPARAM=executable -D CTD_PATH=${CTD_TP_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  # LuciPhorAdapter
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=LuciphorAdapter -DPARAM=executable -D CTD_PATH=${CTD_TP_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  # CometAdapter
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=CometAdapter -DPARAM=comet_executable -D CTD_PATH=${CTD_TP_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  # PercolatorAdapter
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=PercolatorAdapter -DPARAM=percolator_executable -D CTD_PATH=${CTD_TP_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
  DEPENDS create_ctds
)

# remove out_type parameters
foreach(TOOL ${CTD_executables})
  if(${TOOL} IN_LIST THIRDPARTY_ADAPTERS)
    add_custom_command(
      TARGET final_ctds POST_BUILD
      COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=${TOOL} -DPARAM=out_type -D CTD_PATH=${CTD_TP_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
    )
  else()
    add_custom_command(
        TARGET final_ctds POST_BUILD
        COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -DTOOLNAME=${TOOL} -DPARAM=out_type -D CTD_PATH=${CTD_PATH} -P ${SCRIPT_DIRECTORY}remove_parameter_from_ctd.cmake
    )
  endif()
endforeach()

# copy the icons
add_custom_target(
    create_icons
    COMMAND ${CMAKE_COMMAND} -E make_directory ${KNIME_PLUGIN_DIRECTORY}/icons
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/cmake/knime/icons ${KNIME_PLUGIN_DIRECTORY}/icons
    COMMAND ${CMAKE_COMMAND} -E make_directory ${KNIME_TP_PLUGIN_DIRECTORY}/icons
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/icons/category.png ${KNIME_TP_PLUGIN_DIRECTORY}/icons/
    DEPENDS create_knime_folders
)

# TODO in theory thirdparty does not need ALL the mimetypes. Should be ok for now.
# TODO check why we need two files.
add_custom_target(
  prepare_knime_descriptors
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/mimetypes.xml ${CTD_PATH}
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/mime.types ${CTD_PATH}
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/mimetypes.xml ${CTD_TP_PATH}
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/mime.types ${CTD_TP_PATH}
  DEPENDS create_knime_folders final_ctds
)

add_custom_target(
  prepare_knime_payload_binaries
  DEPENDS create_knime_folders
)

# copy the binaries
foreach(TOOL ${CTD_executables})
  set(tool_path ${TOPP_BIN_PATH}/${TOOL}${CMAKE_EXECUTABLE_SUFFIX})
  if(${TOOL} IN_LIST THIRDPARTY_ADAPTERS)
    add_custom_command(
        TARGET  prepare_knime_payload_binaries POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${tool_path} "${TP_PAYLOAD_BIN_PATH}/"
    )
  else()
    add_custom_command(
        TARGET  prepare_knime_payload_binaries POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${tool_path} "${PAYLOAD_BIN_PATH}/"
    )
  endif()
endforeach()

add_custom_target(
    create_payload_share
    # 1st create the directory to make sure that the remove_directory does not fail
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_SHARE_PATH}
    # remove old directory
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PAYLOAD_SHARE_PATH}
    # create new one and fill with the appropriate content
    COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_SHARE_PATH}
    COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SOURCE_PATH=${PROJECT_SOURCE_DIR} -D TARGET_DIRECTORY=${PAYLOAD_SHARE_PATH} -P ${SCRIPT_DIRECTORY}copy_share.cmake
    DEPENDS create_knime_folders
)

add_custom_target(
  prepare_knime_payload_libs
  # we may need the binaries to determine what libraries we need
  DEPENDS prepare_knime_payload_binaries
)

## Kept for reference: we removed the qsqlite dependency and don't need the plugin anymore
## Copy Sqlite plugin and create qt.conf
# TODO in theory we should make that dependent on if Qt was linked dynamically but this is all we support
#  mid-term anyway.
#add_custom_command(
#    TARGET prepare_knime_payload_libs POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E make_directory ${PAYLOAD_LIB_PATH}/plugins/
#    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Qt6::QSQLiteDriverPlugin> ${PAYLOAD_LIB_PATH}/plugins/
#)
# create qt.conf file that specifies plugin dir location
#add_custom_command(
#    TARGET prepare_knime_payload_libs POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E echo "[Paths]" > ${PAYLOAD_BIN_PATH}/qt.conf
#    COMMAND ${CMAKE_COMMAND} -E echo "Plugins = ../${PAYLOAD_LIB_PATH}/plugins" >> ${PAYLOAD_BIN_PATH}/qt.conf
#    COMMAND ${CMAKE_COMMAND} -E echo "" >> ${PAYLOAD_BIN_PATH}/qt.conf
#)

## Assemble common required libraries for win and lnx
## Note that we do not need QtGui libraries since we do not include GUI tools here in KNIME.
## see also the REMOVE_ITEM variable which removes tools with non-obvious dependencies to QtGui
foreach (KNIME_TOOLS_DEPENDENCY OpenMS OpenSwathAlgo)
  add_custom_command(
      TARGET prepare_knime_payload_libs POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${KNIME_TOOLS_DEPENDENCY}> ${PAYLOAD_LIB_PATH}
  )
endforeach()

# assemble the libraries
if (APPLE) ## On APPLE use our script because the executables' install_names need to be changed. Probably can be changed as soon as all
  ## of our dynamically built dependencies build with rpath enabled on brew. Qt recently did the switch for example. This is because if the default install_name of
  ## Qt is /usr/local/qt6/QtCore, then this will be hardcoded in our libOpenMS and tools, even if we use @rpath throughout all of our
  ## cmake build system. See e.g., https://discourse.cmake.org/t/how-to-get-an-lc-rpath-and-rpath-prefix-on-a-dylib-on-macos/5540
  add_custom_command(
    TARGET prepare_knime_payload_libs POST_BUILD
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_dependencies.rb -l ${PAYLOAD_LIB_PATH} -b ${PAYLOAD_BIN_PATH} -f -e "@rpath" -n
    COMMAND ${PROJECT_SOURCE_DIR}/cmake/MacOSX/fix_dependencies.rb -l ${PAYLOAD_LIB_PATH} -b ${TP_PAYLOAD_BIN_PATH} -f -e "@rpath" -n
  ) # -p ${PAYLOAD_LIB_PATH}/plugins not applicable for now
  add_custom_command(
          TARGET prepare_knime_payload_libs POST_BUILD
          COMMAND find "${PAYLOAD_BIN_PATH}" -maxdepth 1 -type f -exec ${CMAKE_STRIP} -S {} "\;"
          COMMAND find "${TP_PAYLOAD_BIN_PATH}" -maxdepth 1 -type f -exec ${CMAKE_STRIP} -S {} "\;"
  )
  add_custom_command(
          TARGET prepare_knime_payload_libs POST_BUILD
          COMMAND find "${PAYLOAD_LIB_PATH}" -type f -name "*.dylib" -exec ${CMAKE_STRIP} -x {} "\;"
          COMMAND find "${PAYLOAD_LIB_PATH}" -type f -name "Qt*" -exec ${CMAKE_STRIP} -x {} "\;"
  )
elseif(WIN32)
  # on Win everything should be linked statically for distribution except Qt
  foreach (KNIME_TOOLS_QT6_DEPENDENCY ${OpenMS_QT_COMPONENTS})
    add_custom_command(
		TARGET prepare_knime_payload_libs POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Qt6::${KNIME_TOOLS_QT6_DEPENDENCY}> ${PAYLOAD_LIB_PATH}
	  )
  endforeach()
  # copying multiple files is possible since CMake 3.5. Last entry is destination. Copy all runtime libs
  # figured out by the CMake InstallRequiredSystemLibraries module
  add_custom_command(
      TARGET prepare_knime_payload_libs POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS} ${PAYLOAD_LIB_PATH}
  )
else()
  foreach (KNIME_DEPENDENCY OpenMS OpenSwathAlgo)
    # copy the dependencies of our libs.
    add_custom_command(
      TARGET prepare_knime_payload_libs POST_BUILD
      COMMAND ${CMAKE_COMMAND} -V -DDEPS="$<TARGET_FILE:${KNIME_DEPENDENCY}>" -DTARGET="${PAYLOAD_LIB_PATH}" -DLOOKUP_DIRS="${OPENMS_CONTRIB_LIBS}/lib\;${QT_INSTALL_BINS}\;${QT_INSTALL_LIBS}" -P ${SCRIPT_DIRECTORY}knime_copy_deps.cmake
      )
  endforeach()
  add_custom_command(
          TARGET prepare_knime_payload_libs POST_BUILD
          COMMAND find "${PAYLOAD_BIN_PATH}" -maxdepth 1 -type f -exec ${CMAKE_STRIP} -s {} "\;"
          COMMAND find "${TP_PAYLOAD_BIN_PATH}" -maxdepth 1 -type f -exec ${CMAKE_STRIP} -s {} "\;"
  )
  add_custom_command(
          TARGET prepare_knime_payload_libs POST_BUILD
          COMMAND find "${PAYLOAD_LIB_PATH}" -type f -name "*.so" -exec ${CMAKE_STRIP} -x {} "\;"
  )
endif()

# handle the binaries.ini
add_custom_target(
  prepare_knime_payload_ini
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D ARCH=${ARCH} -D PLATFORM=${PLATFORM} -D TARGET_DIR=${PAYLOAD_PATH} -D TEMPLATE_FOLDER=${SCRIPT_DIRECTORY} -P ${SCRIPT_DIRECTORY}copy_binaries_ini.cmake
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D ARCH=${ARCH} -D PLATFORM=${PLATFORM} -D TARGET_DIR=${TP_PAYLOAD_PATH} -D TEMPLATE_FOLDER=${SCRIPT_DIRECTORY} -P ${SCRIPT_DIRECTORY}copy_binaries_ini.cmake
  DEPENDS prepare_knime_payload_binaries
)

set(FOLDER_STRUCTURE_MESSAGE "You can clone all third-party binaries from our OpenMS/THIRDPARTY Git submodule/repository but you have to flatten the folder structure such that it is only one level deep with the versions specific for your platform. Do not change the folder names.")

# check if we have valid search engines
if(NOT EXISTS ${SEARCH_ENGINES_DIRECTORY})

  message(WARNING "SEARCH_ENGINES_DIRECTORY not specified or found. Will proceed without shipping third-party executables.
If this is unintended, please specify the path to the search engines to build the KNIME packages and make sure it exists.
${FOLDER_STRUCTURE_MESSAGE}
Then call cmake again with cmake -D SEARCH_ENGINES_DIRECTORY=<Path-To-Checkedout-SE>.")

  # add dummy target
  add_custom_target(
          prepare_knime_payload_searchengines
  )

elseif(NOT EXISTS ${SEARCH_ENGINES_DIRECTORY}/Comet OR NOT EXISTS ${SEARCH_ENGINES_DIRECTORY}/Percolator)
  message(FATAL_ERROR "The given SEARCH_ENGINES_DIRECTORY folder seems to have an invalid layout. ${FOLDER_STRUCTURE_MESSAGE}")
else()
  add_custom_target(
          prepare_knime_payload_searchengines
          COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D SE_PATH=${SEARCH_ENGINES_DIRECTORY} -D TARGET_DIRECTORY=${TP_PAYLOAD_BIN_PATH} -P ${SCRIPT_DIRECTORY}copy_searchengines.cmake
          # We need the folder layout from the bin target
          DEPENDS prepare_knime_payload_binaries
  )
endif()


# the complete payload target
add_custom_target(
  prepare_knime_payload
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D ARCH=${ARCH} -D PLATFORM=${PLATFORM} -D PAYLOAD_FOLDER=${PAYLOAD_PATH} -P ${SCRIPT_DIRECTORY}compress_payload.cmake
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D ARCH=${ARCH} -D PLATFORM=${PLATFORM} -D PAYLOAD_FOLDER=${TP_PAYLOAD_PATH} -P ${SCRIPT_DIRECTORY}compress_payload.cmake
  COMMAND ${CMAKE_COMMAND} -D SCRIPT_DIR=${SCRIPT_DIRECTORY} -D ARCH=${ARCH} -D PLATFORM=${PLATFORM} -D PAYLOAD_FOLDER=${LIB_PATH} -P ${SCRIPT_DIRECTORY}compress_payload.cmake

  DEPENDS prepare_knime_payload_binaries prepare_knime_payload_libs create_payload_share create_icons prepare_knime_payload_ini prepare_knime_payload_searchengines
)

add_custom_target(
  prepare_meta_information
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/LICENSE ${OPENMS_FEATURE_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/COPYRIGHT ${OPENMS_FEATURE_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/DESCRIPTION ${OPENMS_FEATURE_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/LICENSE ${KNIME_PLUGIN_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/COPYRIGHT ${KNIME_PLUGIN_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/DESCRIPTION ${KNIME_PLUGIN_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/LICENSE ${KNIME_LIB_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/COPYRIGHT ${KNIME_LIB_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/DESCRIPTION ${KNIME_LIB_DIRECTORY}/
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/TPLICENSE ${OPENMS_TP_FEATURE_DIRECTORY}/LICENSE
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/TPCOPYRIGHT ${OPENMS_TP_FEATURE_DIRECTORY}/COPYRIGHT
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/TPDESCRIPTION ${OPENMS_TP_FEATURE_DIRECTORY}/DESCRIPTION
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/TPLICENSE ${KNIME_TP_PLUGIN_DIRECTORY}/LICENSE
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/TPCOPYRIGHT ${KNIME_TP_PLUGIN_DIRECTORY}/COPYRIGHT
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/TPDESCRIPTION ${KNIME_TP_PLUGIN_DIRECTORY}/DESCRIPTION
  DEPENDS create_knime_folders
) ## TODO maybe create a separate license file for the Thirdparties on-the-fly by using the actual licenses in the THIRDPARTY submodule

add_custom_target(
  prepare_knime_package
  DEPENDS prepare_meta_information configure_plugin_properties prepare_knime_descriptors prepare_knime_payload
)

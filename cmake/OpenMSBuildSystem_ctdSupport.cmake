# CMake script to generate CTDs for all TOPP tools and UTILS
# Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
# ----------------------------------------------------------
#           OpenMS Mass Spectrometry Framework
# ----------------------------------------------------------
# Maintainer: Stephan Aiche

# create the target directory
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/ctds/)

# create plugin.properties file
configure_file(${PROJECT_SOURCE_DIR}/cmake/knime/plugin.properties.in 
               ${PROJECT_BINARY_DIR}/ctds/plugin.properties)

# copy the icons
file(COPY        ${PROJECT_SOURCE_DIR}/cmake/knime/icons
     DESTINATION ${PROJECT_BINARY_DIR}/ctds/)

# path were the CTDs will be stored
set(CTD_PATH ${PROJECT_BINARY_DIR}/ctds/descriptors)

# path were the executables can be found
set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR})

# list of all tools that can generate CTDs
set(CTD_executables ${TOPP_executables} ${UTILS_executables})

# remove tools that do not produce CTDs
list(REMOVE_ITEM CTD_executables PhosphoScoring OpenMSInfo FuzzyDiff GenericWrapper InspectAdapter MascotAdapter PILISIdentification PILISModelCV PILISModelTrainer PILISSpectraGenerator SvmTheoreticalSpectrumGeneratorTrainer)

# create final target that collects all sub-calls
add_custom_target(
	prepare_knime_package
  COMMAND cmake -E remove_directory ${CTD_PATH}
  COMMAND cmake -E make_directory ${CTD_PATH}
  COMMAND cmake -E copy ${PROJECT_SOURCE_DIR}/cmake/knime/mimetypes.xml ${CTD_PATH}
  DEPENDS TOPP UTILS
)

# call the tools
foreach(TOOL ${CTD_executables})
	add_custom_command(
		TARGET  prepare_knime_package POST_BUILD
		COMMAND ${TOPP_BIN_PATH}/${TOOL} -write_ctd ${CTD_PATH}
    COMMENT "Creating ctd for ${TOOL}"
	)
endforeach()

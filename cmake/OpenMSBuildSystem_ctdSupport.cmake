# CMake script to generate CTDs for all TOPP tools and UTILS
# Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
# ----------------------------------------------------------
#           OpenMS Mass Spectrometry Framework
# ----------------------------------------------------------
# Maintainer: Stephan Aiche

# Path were the CTDs will be stored
set(CTD_PATH ${PROJECT_BINARY_DIR}/ctds/descriptors)
file(MAKE_DIRECTORY ${CTD_PATH})

# Path were the executables can be found
set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR})

# list of all tools that can generate CTDs
set(CTD_executables ${TOPP_executables} ${UTILS_executables})
# remove tools that do not produce CTDs
list(REMOVE_ITEM CTD_executables PhosphoScoring OpenMSInfo)

set(CTD_TARGETs "")

# call the tools
FOREACH(TOOL ${CTD_executables})
	add_custom_target(
		${TOOL}_CTD
		COMMAND
		${TOPP_BIN_PATH}/${TOOL} -write_ctd ${CTD_PATH}
	)
	list(APPEND CTD_TARGETs ${TOOL}_CTD)
ENDFOREACH()

# create final target that collects all sub-calls
add_custom_target(
	build_ctds
	DEPENDS
	${CTD_TARGETs}
)

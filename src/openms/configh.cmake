# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Chris Bielow $
# $Authors: Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

## define some source directories
set(CF_OPENMS_DATA_PATH ${OPENMS_HOST_DIRECTORY}/share/OpenMS CACHE INTERNAL "Path to the shared documents of OpenMS.")
set(CF_OPENMS_DOC_PATH ${OPENMS_HOST_DIRECTORY}/doc CACHE INTERNAL "Path to the documentation of OpenMS.")
## and the corresponding ones when installed (careful, you have to rebuild if you change -DCMAKE_PREFIX_PATH). Also, does not work after deployment.
set(CF_OPENMS_INSTALL_DATA_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_SHARE_DIR} CACHE INTERNAL "Path to the installed shared documents of OpenMS.")
set(CF_OPENMS_INSTALL_DOC_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_DOC_DIR} CACHE INTERNAL "Path to the installed documentation of OpenMS." )

#------------------------------------------------------------------------------
# At this point make a summary of where data and doc will be located:
message(STATUS "Info: CF_OPENMS_DATA_PATH: ${CF_OPENMS_DATA_PATH}")
message(STATUS "Info: CF_OPENMS_DOC_PATH: ${CF_OPENMS_DOC_PATH}")


## check for Microsoft Visual Studio compiler
if (MSVC)
	set(OPENMS_COMPILER_MSVC "1" CACHE INTERNAL "Do we use Microsoft Compiler?")
endif()
## check for G++
if (CMAKE_COMPILER_IS_GNUCXX)
	set(OPENMS_COMPILER_GXX "1" CACHE INTERNAL "Do we use G++ Compiler?")
endif()

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(OPENMS_BIG_ENDIAN)

## check 32/64 bit architecture (defined above!)
if (NOT DEFINED OPENMS_64BIT_ARCHITECTURE)
	message(FATAL_ERROR "Cmake script was re-ordered and is now invalid! Please make sure that OPENMS_64BIT_ARCHITECTURE is defined when config.h.in is configured!")
endif()

## conditionally include //@dot commands in doxygen using using #ifdef OPENMS_HASDOXYGENDOT
if (DOXYGEN_HAVE_DOT)
  set(CF_OPENMS_HASDOXYGENDOT 1)
else()
  set(CF_OPENMS_HASDOXYGENDOT 0)
endif()

#------------------------------------------------------------------------------
## Check if various system headers exist
include(CheckIncludeFileCXX)

CHECK_INCLUDE_FILE_CXX("unistd.h" OPENMS_HAS_UNISTD_H)
CHECK_INCLUDE_FILE_CXX("process.h" OPENMS_HAS_PROCESS_H)

CHECK_INCLUDE_FILE_CXX("time.h" OPENMS_HAS_TIME_H)
CHECK_INCLUDE_FILE_CXX("sys/times.h" OPENMS_HAS_SYS_TIMES_H)
CHECK_INCLUDE_FILE_CXX("sys/time.h"  OPENMS_HAS_SYS_TIME_H)
CHECK_INCLUDE_FILE_CXX("sys/resource.h"  OPENMS_HAS_SYS_RESOURCE_H)

#------------------------------------------------------------------------------
# check if certain c++ functions exist
include(CheckFunctionExists)
## in MinGW we have the signal.h header, but no kill() as in Linux, so we need to check for the kill() function
CHECK_FUNCTION_EXISTS("kill" OPENMS_HAS_KILL)
CHECK_FUNCTION_EXISTS("sysconf" OPENMS_HAS_SYSCONF)


#------------------------------------------------------------------------------
# Create the config.h
# replace any variables in config.h.in with current values
set (CONFIGURED_CONFIG_H ${PROJECT_BINARY_DIR}/include/OpenMS/config.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/config.h.in ${CONFIGURED_CONFIG_H})

#------------------------------------------------------------------------------
# Create openms_package_version.h
# replace any variables in openms_package_version.h.in with current values
set (CONFIGURED_OPENMS_PACKAGE_VERSION_H ${PROJECT_BINARY_DIR}/include/OpenMS/openms_package_version.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/openms_package_version.h.in ${CONFIGURED_OPENMS_PACKAGE_VERSION_H})

#------------------------------------------------------------------------------
# create paths header
set(CONFIGURED_OPENMS_DATA_PATH_H ${PROJECT_BINARY_DIR}/include/OpenMS/openms_data_path.h)
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/openms_data_path.h.in ${CONFIGURED_OPENMS_DATA_PATH_H})

#------------------------------------------------------------------------------
# Generate build_config_$config.h at configure time
# Modify the used build_config.h at build time
set (CONFIGURED_BUILD_CONFIG_H ${PROJECT_BINARY_DIR}/include/OpenMS/build_config_$<CONFIG>.h)
set (CONFIGURED_BUILD_CONFIG_CURRENT_H ${PROJECT_BINARY_DIR}/include/OpenMS/build_config.h)

file (GENERATE
			OUTPUT ${CONFIGURED_BUILD_CONFIG_H}
			CONTENT "
#ifndef OPENMS_BUILD_CONFIG_H
#define OPENMS_BUILD_CONFIG_H

#define OPENMS_BUILD_TYPE \"$<CONFIG>\"

#endif // OPENMS_BUILD_CONFIG_H
"
			)

add_custom_command (
		COMMAND ${CMAKE_COMMAND} "-E" "copy_if_different" "${CONFIGURED_BUILD_CONFIG_H}" "${CONFIGURED_BUILD_CONFIG_CURRENT_H}"
		VERBATIM
		DEPENDS  "${CONFIGURED_BUILD_CONFIG_H}"
		OUTPUT   "${CONFIGURED_BUILD_CONFIG_CURRENT_H}"
		COMMENT  "creating build_config.h file ({event: PRE_BUILD}, {filename: build_config.h})"
)

#------------------------------------------------------------------------------
# export a list of all configured headers
set(OpenMS_configured_headers "${CONFIGURED_BUILD_CONFIG_CURRENT_H};${CONFIGURED_CONFIG_H};${CONFIGURED_OPENMS_PACKAGE_VERSION_H};${CONFIGURED_OPENMS_DATA_PATH_H}")
set_property(SOURCE ${CONFIGURED_BUILD_CONFIG_CURRENT_H} PROPERTY SKIP_AUTOMOC ON)

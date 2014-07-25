# - Try to find SEQAN headers
# Once done, this will define
#
#  SEQAN_FOUND - system has SEQAN headers installed (no libraries required)
#  SEQAN_INCLUDE_DIRS - the SEQAN include directories

# Include dir
find_path(SEQAN_INCLUDE_DIR seqan/version.h)
if (SEQAN_INCLUDE_DIR)
	file(READ "${SEQAN_INCLUDE_DIR}/seqan/version.h" _SEQAN_VERSION_H_CONTENTS)
	string(REGEX REPLACE ".*#define SEQAN_VERSION_MAJOR ([0-9]+).*" "\\1" SEQAN_VERSION_MAJOR "${_SEQAN_VERSION_H_CONTENTS}")
	string(REGEX REPLACE ".*#define SEQAN_VERSION_MINOR ([0-9]+).*" "\\1" SEQAN_VERSION_MINOR "${_SEQAN_VERSION_H_CONTENTS}")
	string(REGEX REPLACE ".*#define SEQAN_VERSION_PATCH ([0-9]+).*" "\\1" SEQAN_VERSION_PATCH "${_SEQAN_VERSION_H_CONTENTS}")
  set(SEQAN_VERSION_STRING "${SEQAN_VERSION_MAJOR}.${SEQAN_VERSION_MINOR}.${SEQAN_VERSION_PATCH}")
endif()

set(SEQAN_INCLUDE_DIRS ${SEQAN_INCLUDE_DIR})

#------------------------------------------------------------------------------
# handle standard args
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SEQAN
                                  REQUIRED_VARS SEQAN_INCLUDE_DIRS
                                  VERSION_VAR SEQAN_VERSION_STRING)

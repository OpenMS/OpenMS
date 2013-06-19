# - Try to find SEQAN headers
# Once done, this will define
#
#  SEQAN_FOUND - system has SEQAN headers installed (no libraries required)
#  SEQAN_INCLUDE_DIRS - the SEQAN include directories

include(LibFindMacros)

# Include dir
find_path(SEQAN_INCLUDE_DIR seqan/version.h)
if (NOT SEQAN_INCLUDE_DIR)
	MESSAGE(STATUS "Could not find SEQAN installation.")
	set(SEQAN_FOUND 0)
else()
	set(SEQAN_FOUND 1)
	file(READ "${SEQAN_INCLUDE_DIR}/seqan/version.h" _SEQAN_VERSION_H_CONTENTS)
	string(REGEX REPLACE ".*#define SEQAN_VERSION_MAJOR ([0-9]+).*" "\\1" SEQAN_VERSION_MAJOR "${_SEQAN_VERSION_H_CONTENTS}")
	string(REGEX REPLACE ".*#define SEQAN_VERSION_MINOR ([0-9]+).*" "\\1" SEQAN_VERSION_MINOR "${_SEQAN_VERSION_H_CONTENTS}")
	string(REGEX REPLACE ".*#define SEQAN_VERSION_PATCH ([0-9]+).*" "\\1" SEQAN_VERSION_PATCH "${_SEQAN_VERSION_H_CONTENTS}")

	# Set the include dir variables and the libraries and let libfind_process do the rest.
	# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
	set(SEQAN_PROCESS_INCLUDES SEQAN_INCLUDE_DIR SEQAN_INCLUDE_DIRS)
	set(SEQAN_PROCESS_LIBS SEQAN_LIBRARY SEQAN_LIBRARIES)
endif()




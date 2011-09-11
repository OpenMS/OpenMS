# - Try to find SEQAN headers
# Once done, this will define
#
#  SEQAN_FOUND - system has SEQAN headers installed (no libraries required)
#  SEQAN_INCLUDE_DIRS - the SEQAN include directories

include(LibFindMacros)

# Include dir
find_path(SEQAN_INCLUDE_DIR seqan.h)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(SEQAN_PROCESS_INCLUDES SEQAN_INCLUDE_DIR SEQAN_INCLUDE_DIRS)
set(SEQAN_PROCESS_LIBS SEQAN_LIBRARY SEQAN_LIBRARIES)

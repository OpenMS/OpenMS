# - Try to find libssh using pkg-config
# Once done, this will define
#
#  LIBSSH2_FOUND - system has libssh2
#  LIBSSH2_INCLUDE_DIRS - the libssh2 include directories
#  LIBSSH2_LIBRARIES - link these to use libssh2

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(LIBSSH2_PKGCONF libssh2)

# Include dir
find_path(LIBSSH2_INCLUDE_DIR
  NAMES libssh2.h
  PATHS LIBSSH2_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(LIBSSH2_LIBRARY
  NAMES ssh2
  PATHS ${LIBSSH2_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(LIBSSH2_PROCESS_INCLUDES LIBSSH2_INCLUDE_DIR LIBSSH2_INCLUDE_DIRS)
set(LIBSSH2_PROCESS_LIBS LIBSSH2_LIBRARY LIBSSH2_LIBRARIES)
libfind_process(LIBSSH2)

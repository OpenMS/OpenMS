# - Try to find libssh using pkg-config
# Once done, this will define
#
#  LIBSSH_FOUND - system has libssh
#  LIBSSH_INCLUDE_DIRS - the libssh include directories
#  LIBSSH_LIBRARIES - link these to use libssh

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(LIBSSH_PKGCONF libssh)

# Include dir
find_path(LIBSSH_INCLUDE_DIR
  NAMES libssh.h
  PATHS LIBSSH_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(LIBSSH_LIBRARY
  NAMES ssh
  PATHS ${LIBSSH_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(LIBSSH_PROCESS_INCLUDES LIBSSH_INCLUDE_DIR LIBSSH_INCLUDE_DIRS)
set(LIBSSH_PROCESS_LIBS LIBSSH_LIBRARY LIBSSH_LIBRARIES)
libfind_process(LIBSSH)

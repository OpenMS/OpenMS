# - Try to find BZip2
# Once done this will define
#
#  BZIP2_FOUND - system has BZip2
#  BZIP2_INCLUDE_DIR - the BZip2 include directory
#  BZIP2_LIBRARIES - Link these to use BZip2
#  BZIP2_NEED_PREFIX - this is set if the functions are prefixed with BZ2_
#  BZIP2_VERSION_STRING - the version of BZip2 found (since CMake 2.8.8)

#=============================================================================
# Copyright 2006-2012 Kitware, Inc.
# Copyright 2006 Alexander Neundorf <neundorf@kde.org>
# Copyright 2012 Rolf Eike Beer <eike@sf-mail.de>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

set(_BZIP2_PATHS PATHS
  "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\Bzip2;InstallPath]"
  )

find_path(KISSFFT_INCLUDE_DIR kiss_fft.h ${_BZIP2_PATHS} PATH_SUFFIXES include)

if (NOT KISSFFT_LIBRARIES)
    find_library(KISSFFT_LIBRARY_RELEASE NAMES libkfft ${_KISFFT_PATHS} PATH_SUFFIXES lib)
    find_library(KISSFFT_LIBRARY_DEBUG NAMES libkfftd ${_KISSFFT_PATHS} PATH_SUFFIXES lib)

    include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)
    SELECT_LIBRARY_CONFIGURATIONS(KISSFFT)
endif ()

if (KISSFFT_INCLUDE_DIR AND EXISTS "${KISSFFT_INCLUDE_DIR}/kiss_fft.h")
    file(STRINGS "${KISSFFT_INCLUDE_DIR}/kiss_fft.h" KISS_FFT_H REGEX "kissfft/libkfft version [0-9]+\\.[^ ]+ of [0-9]+ ")
    string(REGEX REPLACE ".* kissfft/libkfft version ([0-9]+\\.[^ ]+) of [0-9]+ .*" "\\1" KISSFFT_VERSION_STRING "${KISS_FFT_H}")
endif ()

# handle the QUIETLY and REQUIRED arguments and set BZip2_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(KissFFT
                                  REQUIRED_VARS KISSFFT_LIBRARIES KISSFFT_INCLUDE_DIR
                                  VERSION_VAR KISSFFT_VERSION_STRING)

if (KISSFFT_FOUND)
   include(CheckLibraryExists)
#   CHECK_LIBRARY_EXISTS("${KISSFFT_LIBRARIES}" BZ2_bzCompressInit "" BZIP2_NEED_PREFIX)
endif ()

mark_as_advanced(KISSFFT_INCLUDE_DIR)

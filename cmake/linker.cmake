# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Choice of linker (ELF platforms with GCC or Clang)
#
# libOpenMS and ~900 executables (TOPP tools, class tests) are linked in a full
# build, and every change to libOpenMS relinks all of them. The GNU BFD linker
# takes ~6x longer for this than lld or mold, so a fast linker is used when one
# is available.
#
# OPENMS_LINKER:
#   auto    (default) lld, or mold for GCC with ENABLE_IPO (lld cannot read GCC
#           LTO objects); falls back to the toolchain default if neither works
#   lld | mold | gold | bfd   use exactly this linker (error if it does not work)
#   default the toolchain's default linker (no -fuse-ld)
#
# Nothing is changed when CMAKE_LINKER_TYPE (CMake >= 3.29) or a -fuse-ld= in the
# linker flags already picks a linker.
#------------------------------------------------------------------------------
set(OPENMS_LINKER "auto" CACHE STRING "Linker to use: auto, lld, mold, gold, bfd or default (see cmake/linker.cmake)")
set_property(CACHE OPENMS_LINKER PROPERTY STRINGS auto lld mold gold bfd default)

set(_openms_linker_flags "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_MODULE_LINKER_FLAGS}")
if(OPENMS_LINKER STREQUAL "default")
  message(STATUS "Linker: toolchain default (OPENMS_LINKER=default)")
elseif(APPLE OR MSVC OR NOT (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
  message(STATUS "Linker: toolchain default (OPENMS_LINKER only applies to ELF platforms with GCC/Clang)")
elseif(CMAKE_LINKER_TYPE)
  message(STATUS "Linker: ${CMAKE_LINKER_TYPE} (from CMAKE_LINKER_TYPE)")
elseif(_openms_linker_flags MATCHES "-fuse-ld=")
  message(STATUS "Linker: chosen by -fuse-ld= in the CMAKE_*_LINKER_FLAGS")
else()
  include(CheckLinkerFlag)
  if(OPENMS_LINKER STREQUAL "auto")
    if(ENABLE_IPO AND CMAKE_COMPILER_IS_GNUCXX)
      set(_openms_linker_candidates mold)
    else()
      set(_openms_linker_candidates lld mold)
    endif()
  else()
    set(_openms_linker_candidates ${OPENMS_LINKER})
  endif()

  set(_openms_linker_chosen "")
  foreach(_ld IN LISTS _openms_linker_candidates)
    # the result is cached per linker; a changed toolchain needs a fresh cache
    check_linker_flag(CXX "-fuse-ld=${_ld}" OPENMS_LINKER_SUPPORTS_${_ld})
    if(OPENMS_LINKER_SUPPORTS_${_ld})
      set(_openms_linker_chosen ${_ld})
      break()
    endif()
  endforeach()

  if(_openms_linker_chosen)
    add_link_options("-fuse-ld=${_openms_linker_chosen}")
    message(STATUS "Linker: ${_openms_linker_chosen} (OPENMS_LINKER=${OPENMS_LINKER})")
  elseif(OPENMS_LINKER STREQUAL "auto")
    message(STATUS "Linker: toolchain default (no lld/mold found; installing one speeds up linking considerably)")
  else()
    message(FATAL_ERROR "OPENMS_LINKER=${OPENMS_LINKER}, but the compiler cannot link with -fuse-ld=${OPENMS_LINKER}. "
                        "Install it, or use OPENMS_LINKER=auto or default.")
  endif()
endif()
unset(_openms_linker_flags)

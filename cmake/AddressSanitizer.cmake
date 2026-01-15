# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Hannes Roest $
# $Authors: Hannes Roest, Timo Sachsenberg $
# --------------------------------------------------------------------------


#------------------------------------------------------------------------------
# This cmake file enables the AddressSanitizer and UndefinedBehaviorSanitizer
# see http://clang.llvm.org/docs/AddressSanitizer.html and https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html
#     http://en.wikipedia.org/wiki/AddressSanitizer

# Do the check for support in the beginning and NOT FOR EVERY TARGET
#------------------------------------------------------------------------------
if(ADDRESS_SANITIZER)
  if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    # add compiler flag
    if (MSVC AND MSVC_VERSION LESS 1920)
        message(FATAL_ERROR "AddressSanitizer is not supported for MSVC versions < 2019.")
    endif()
  else()
    message(FATAL_ERROR "AddressSanitizer is supported for OpenMS debug mode only.")
  endif()
  message(STATUS "AddressSanitizer enabled. Adding flags to every target.")
endif()


function(add_asan_to_target TARGET_NAME_ARG)
  if(ADDRESS_SANITIZER)
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
      # add compiler flag
      if (MSVC)
        # add AddressSanitizer for compiler and linker
        target_compile_options("${TARGET_NAME_ARG}" 
          PUBLIC 
            /fsanitize=address
            )
        target_link_options("${TARGET_NAME_ARG}" 
          PUBLIC 
            /fsanitize=address
            )
      else()
        # add AddressSanitizer for compiler and linker
        target_compile_options("${TARGET_NAME_ARG}" 
          PUBLIC 
            -fsanitize=address,undefined
            -fno-sanitize-recover=all
            -fno-sanitize=vptr
            -fno-omit-frame-pointer)
        target_link_options("${TARGET_NAME_ARG}" 
          PUBLIC 
            -fsanitize=address,undefined
            -fno-sanitize-recover=all
            -fno-sanitize=vptr)
      endif()
    endif()
  endif()
endfunction()

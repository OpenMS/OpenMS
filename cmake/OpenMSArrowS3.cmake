# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause

# ArrowConfig.cmake imports these capabilities from the installed ArrowOptions.cmake.
# They describe that installation; setting them here cannot add S3 support to it.
function(openms_check_arrow_s3)
  if(NOT WITH_S3)
    return()
  endif()
  if(NOT ARROW_S3)
    message(FATAL_ERROR
      "WITH_S3=ON requires Apache Arrow built with ARROW_S3=ON and "
      "ARROW_S3_MODULE=OFF. Install an Arrow package with built-in S3 support "
      "or configure OpenMS with -DWITH_S3=OFF to disable S3 XML loading. "
      "OpenMS does not find or build a separate AWS SDK.")
  endif()
  if(ARROW_S3_MODULE)
    message(FATAL_ERROR
      "WITH_S3=ON requires S3 built into Apache Arrow (ARROW_S3=ON, "
      "ARROW_S3_MODULE=OFF). Separate dynamic Arrow S3 modules are not supported "
      "by OpenMS in this initial implementation. "
      "Rebuild Arrow with these options or configure OpenMS with -DWITH_S3=OFF.")
  endif()
  message(STATUS "S3 XML loading enabled through Arrow's built-in S3 support")
endfunction()

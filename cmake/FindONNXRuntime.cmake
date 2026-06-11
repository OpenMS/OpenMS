#  Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
#  SPDX-License-Identifier: BSD-3-Clause
#
#  --------------------------------------------------------------------------
#  $Maintainer: Timo Sachsenberg $
#  $Authors: Satyam Yadav $
#  --------------------------------------------------------------------------

# FindONNXRuntime.cmake
#
# Finds the ONNX Runtime library
#
# This will define the following variables:
#   ONNXRuntime_FOUND        - True if the system has the ONNX Runtime library
#   ONNXRuntime_INCLUDE_DIRS - The ONNX Runtime include directories
#   ONNXRuntime_LIBRARIES    - The ONNX Runtime libraries
#
# And the following imported target:
#   ONNXRuntime::ONNXRuntime

find_path(ONNXRuntime_INCLUDE_DIR
  NAMES onnxruntime_cxx_api.h
  PATH_SUFFIXES include onnxruntime/core/session
)

find_library(ONNXRuntime_LIBRARY
  NAMES onnxruntime
  PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ONNXRuntime
  REQUIRED_VARS ONNXRuntime_LIBRARY ONNXRuntime_INCLUDE_DIR
)

if(ONNXRuntime_FOUND)
  set(ONNXRuntime_INCLUDE_DIRS ${ONNXRuntime_INCLUDE_DIR})
  set(ONNXRuntime_LIBRARIES ${ONNXRuntime_LIBRARY})

  if(NOT TARGET ONNXRuntime::ONNXRuntime)
    add_library(ONNXRuntime::ONNXRuntime UNKNOWN IMPORTED)
    set_target_properties(ONNXRuntime::ONNXRuntime PROPERTIES
      IMPORTED_LOCATION "${ONNXRuntime_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ONNXRuntime_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(ONNXRuntime_INCLUDE_DIR ONNXRuntime_LIBRARY)
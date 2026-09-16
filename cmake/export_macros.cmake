# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


include(CMakePackageConfigHelpers)

# a collection of wrapper for export functions that allows easier usage
# througout the OpenMS build system
set(_OPENMS_EXPORT_FILE "OpenMSTargets.cmake")

# clear list before we refill it
set(_OPENMS_EXPORT_TARGETS "" CACHE INTERNAL "List of targets that will be exported.")

macro(openms_register_export_target target_name)
  set(_OPENMS_EXPORT_TARGETS ${_OPENMS_EXPORT_TARGETS} ${target_name}
    CACHE INTERNAL "List of targets that will be exported.")
endmacro()

macro(openms_export_targets )

  # configure OpenMSConfig.cmake
  configure_package_config_file(
    "${OPENMS_HOST_DIRECTORY}/cmake/OpenMSConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/OpenMSConfig.cmake"
    INSTALL_DESTINATION ${INSTALL_CMAKE_DIR}
    PATH_VARS INSTALL_SHARE_DIR INSTALL_LIB_DIR INSTALL_DOC_DIR INSTALL_BIN_DIR
  )

  # write OpenMSConfigVersion.cmake
  write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/OpenMSConfigVersion.cmake"
    VERSION ${OPENMS_PACKAGE_VERSION}
    COMPATIBILITY SameMinorVersion
  )

  # create the corresponding target file for the build tree, with the same
  # OpenMS:: namespace as the installed export (install_export_targets), so a
  # project configured against the build tree sees the same target names
  export(TARGETS ${_OPENMS_EXPORT_TARGETS}
         NAMESPACE OpenMS::
         FILE ${OPENMS_HOST_BINARY_DIRECTORY}/${_OPENMS_EXPORT_FILE})

  # install the generated config file
  install_file(${PROJECT_BINARY_DIR}/OpenMSConfig.cmake
               ${INSTALL_CMAKE_DIR}
               cmake)

  # .. and ConfigVersion.cmake
  install_file(${PROJECT_BINARY_DIR}/OpenMSConfigVersion.cmake
               ${INSTALL_CMAKE_DIR}
               cmake)

  # No export(PACKAGE OpenMS): with cmake_minimum_required(VERSION 3.24) policy
  # CMP0090 makes it a no-op, and consumers select an installation explicitly
  # through OpenMS_DIR or CMAKE_PREFIX_PATH rather than via the user package registry.
endmacro()

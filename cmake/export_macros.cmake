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
#
# The export sets (core, CLI, GUI; see cmake/install_macros.cmake) each get one
# target file, <set>.cmake, in the build tree and in the installation.
# OpenMSConfig.cmake includes the core file and then the files of the layers
# that are present.
set(_OPENMS_EXPORT_FILE "${OPENMS_EXPORT_SET}.cmake")
set(_OPENMS_CLI_EXPORT_FILE "${OPENMS_CLI_EXPORT_SET}.cmake")
set(_OPENMS_GUI_EXPORT_FILE "${OPENMS_GUI_EXPORT_SET}.cmake")

# clear the lists before we refill them (one list of targets per export set)
foreach(_export_set IN LISTS OPENMS_EXPORT_SETS)
  set(_OPENMS_EXPORT_TARGETS_${_export_set} "" CACHE INTERNAL "Targets exported in ${_export_set}.")
endforeach()

# openms_register_export_target(<target> [<export set>])
# Registers a target for the build-tree export of its export set (the core set
# OpenMSTargets when omitted). The install-tree export is registered by
# install_library().
macro(openms_register_export_target target_name)
  set(_export_set ${ARGN})
  if("${_export_set}" STREQUAL "")
    set(_export_set ${OPENMS_EXPORT_SET})
  endif()
  if(NOT "${_export_set}" IN_LIST OPENMS_EXPORT_SETS)
    message(FATAL_ERROR "openms_register_export_target(${target_name}): unknown export set '${_export_set}' (known: ${OPENMS_EXPORT_SETS})")
  endif()
  set(_OPENMS_EXPORT_TARGETS_${_export_set} ${_OPENMS_EXPORT_TARGETS_${_export_set}} ${target_name}
    CACHE INTERNAL "Targets exported in ${_export_set}.")
  unset(_export_set)
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

  # create the corresponding target files for the build tree, one per export
  # set that has targets, with the same OpenMS:: namespace as the installed
  # export (install_export_targets), so a project configured against the build
  # tree sees the same target names. A target of one set refers to targets of
  # another set through that set's file, which CMake resolves across the
  # export() calls of this project.
  foreach(_export_set IN LISTS OPENMS_EXPORT_SETS)
    if(_OPENMS_EXPORT_TARGETS_${_export_set})
      export(TARGETS ${_OPENMS_EXPORT_TARGETS_${_export_set}}
             NAMESPACE OpenMS::
             FILE ${OPENMS_HOST_BINARY_DIRECTORY}/${_export_set}.cmake)
    endif()
  endforeach()

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

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

if(TYPE STREQUAL "LIBPLUGIN")
  set(ID "de.openms.lib")
  set(NAME "OpenMSLibraries")
  # create plugin.properties file
  configure_file(${SOURCE_PATH}/cmake/knime/libplugin.properties.in
                 ${TARGET_PATH}/plugin.properties)
elseif(TYPE STREQUAL "TOPPPLUGIN")
  set(ID "de.openms")
  set(NAME "OpenMS")
  # create plugin.properties file
  configure_file(${SOURCE_PATH}/cmake/knime/plugin.properties.in
                 ${TARGET_PATH}/plugin.properties)
elseif(TYPE STREQUAL "TPPLUGIN")
  set(ID "de.openms.thirdparty")
  set(NAME "OpenMSThirdparty")
  # create plugin.properties file
  configure_file(${SOURCE_PATH}/cmake/knime/plugin.properties.in
                 ${TARGET_PATH}/plugin.properties)
elseif(TYPE STREQUAL "FEATURE")
  set(ID "de.openms.feature")
  set(NAME "OpenMS")
  set(CATEGORY "openms")
  # create plugin.properties file
  configure_file(${SOURCE_PATH}/cmake/knime/feature.properties.in
                 ${TARGET_PATH}/feature.properties)
elseif(TYPE STREQUAL "TPFEATURE")
  set(ID "de.openms.thirdparty.feature")
  set(NAME "OpenMSThirdparty")
  set(CATEGORY "openmsthirdparty")
  # create plugin.properties file
  configure_file(${SOURCE_PATH}/cmake/knime/feature.properties.in
                 ${TARGET_PATH}/feature.properties)
elseif(TYPE STREQUAL "UPDATESITE")
  # create plugin.properties file
  configure_file(${SOURCE_PATH}/cmake/knime/updatesite.properties.in
                 ${TARGET_PATH}/updatesite.properties)
endif()

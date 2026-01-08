# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


project("OpenMS-evergreen")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# export include directory
add_library(Evergreen INTERFACE)
target_include_directories(Evergreen SYSTEM INTERFACE 
                            "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>"
                            "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty/evergreen>")

# Export it under OpenMSTargets.cmake instead of a separate EvergreenTargets file
#install(EXPORT EvergreenTargets DESTINATION lib/cmake/evergreen)
#install(TARGETS Evergreen EXPORT EvergreenTargets PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty/evergreen COMPONENT library)

install_library(Evergreen)
openms_register_export_target(Evergreen)

get_filename_component(EVERGREEN_HEADERS "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(EVERGREEN_FOLDER "${EVERGREEN_HEADERS}" PATH)
file(APPEND ${TPLICENSEFILE} "\nevergreen:\n")
cat(${EVERGREEN_FOLDER}/LICENSE ${TPLICENSEFILE})

# Install headers
install(DIRECTORY ${EVERGREEN_HEADERS}/
  DESTINATION ${INSTALL_INCLUDE_DIR}/OpenMS/thirdparty/evergreen
  COMPONENT thirdparty_headers
  FILE_PERMISSIONS      OWNER_WRITE OWNER_READ
                        GROUP_READ
                        WORLD_READ
  DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_EXECUTE GROUP_READ
                        WORLD_EXECUTE WORLD_READ
  REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
  REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
  PATTERN "demos/*" EXCLUDE
  )

# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: $
# $Authors: $
# --------------------------------------------------------------------------

# This script generates an ApplicationsComponent.plist file containing
# an array of all installed application bundles for the Applications CPack component.
# This is used by the CPack productbuild generator on macOS.

# The plist file will be generated during CMake configure and referenced
# in cpack_add_component(Applications ...) via the PLIST argument.

# Get the list of GUI executables
include(${PROJECT_SOURCE_DIR}/src/openms_gui/source/VISUAL/APPLICATIONS/GUITOOLS/executables.cmake)

# Generate the plist file in the build directory
set(APPLICATIONS_COMPONENT_PLIST "${CMAKE_BINARY_DIR}/ApplicationsComponent.plist")

# Start the plist XML
set(PLIST_CONTENT "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
set(PLIST_CONTENT "${PLIST_CONTENT}<!DOCTYPE plist PUBLIC \"-//Apple//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">\n")
set(PLIST_CONTENT "${PLIST_CONTENT}<plist version=\"1.0\">\n")
set(PLIST_CONTENT "${PLIST_CONTENT}<array>\n")

# Add an entry for each GUI application bundle
foreach(app_name ${GUI_executables})
  set(PLIST_CONTENT "${PLIST_CONTENT}    <dict>\n")
  set(PLIST_CONTENT "${PLIST_CONTENT}        <key>RootRelativeBundlePath</key>\n")
  set(PLIST_CONTENT "${PLIST_CONTENT}        <string>${app_name}.app</string>\n")
  set(PLIST_CONTENT "${PLIST_CONTENT}        <key>BundleIsRelocatable</key>\n")
  set(PLIST_CONTENT "${PLIST_CONTENT}        <true/>\n")
  set(PLIST_CONTENT "${PLIST_CONTENT}    </dict>\n")
endforeach()

# Close the array and plist
set(PLIST_CONTENT "${PLIST_CONTENT}</array>\n")
set(PLIST_CONTENT "${PLIST_CONTENT}</plist>\n")

# Write the plist file
file(WRITE "${APPLICATIONS_COMPONENT_PLIST}" "${PLIST_CONTENT}")
message(STATUS "Generated ApplicationsComponent.plist with bundles: ${GUI_executables}")

# Export the path so it can be used in package_components.cmake
set(APPLICATIONS_COMPONENT_PLIST "${APPLICATIONS_COMPONENT_PLIST}" PARENT_SCOPE)

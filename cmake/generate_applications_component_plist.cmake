# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------

## Writes the component property list that pkgbuild gets for the Applications component
## (PLIST of cpack_add_component in cmake/package_components.cmake) and sets
## APPLICATIONS_COMPONENT_PLIST to it. Included by cmake/package_mac_productbuild.cmake,
## after CPACK_PACKAGING_INSTALL_PREFIX is set.
##
## pkgbuild marks every app bundle relocatable by default. The installer then looks for a
## bundle with the same identifier anywhere on the disk and updates that copy instead of
## installing into the package's folder. Our identifiers (de.openms.<name>) are the same in
## every version, so the apps of a new OpenMS ended up in the folder of an older one, or in
## a build tree. pkgbuild reads BundleIsRelocatable only from this plist, not from the
## Info.plist of a bundle. The other keys keep the values pkgbuild uses by default.

## add_mac_app_bundle() (src/openms_gui/add_mac_bundle.cmake) records every bundle it
## installs; there are none without WITH_GUI.
get_property(_openms_app_bundles GLOBAL PROPERTY OPENMS_APP_BUNDLES)
if(NOT _openms_app_bundles)
  return()
endif()

## The path is relative to the root pkgbuild packages, which holds the install prefix:
## Applications/OpenMS-<version>/TOPPView.app
string(REGEX REPLACE "^/" "" _openms_bundle_dir "${CPACK_PACKAGING_INSTALL_PREFIX}")
## Packages of other branches carry the branch name in the version (CMakeLists.txt), and a
## branch name may contain characters that XML reserves; & first.
string(REPLACE "&" "&amp;" _openms_bundle_dir "${_openms_bundle_dir}")
string(REPLACE "<" "&lt;" _openms_bundle_dir "${_openms_bundle_dir}")
string(REPLACE ">" "&gt;" _openms_bundle_dir "${_openms_bundle_dir}")

set(_openms_plist_entries "")
foreach(_openms_app_bundle IN LISTS _openms_app_bundles)
  string(APPEND _openms_plist_entries
    "  <dict>\n"
    "    <key>BundleHasStrictIdentifier</key>\n"
    "    <true/>\n"
    "    <key>BundleIsRelocatable</key>\n"
    "    <false/>\n"
    "    <key>BundleIsVersionChecked</key>\n"
    "    <true/>\n"
    "    <key>BundleOverwriteAction</key>\n"
    "    <string>upgrade</string>\n"
    "    <key>RootRelativeBundlePath</key>\n"
    "    <string>${_openms_bundle_dir}/${_openms_app_bundle}.app</string>\n"
    "  </dict>\n")
endforeach()

set(APPLICATIONS_COMPONENT_PLIST "${CMAKE_BINARY_DIR}/ApplicationsComponent.plist")
file(WRITE "${APPLICATIONS_COMPONENT_PLIST}"
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
  "<!DOCTYPE plist PUBLIC \"-//Apple//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">\n"
  "<plist version=\"1.0\">\n"
  "<array>\n"
  "${_openms_plist_entries}"
  "</array>\n"
  "</plist>\n")
list(JOIN _openms_app_bundles ", " _openms_app_bundles_text)
message(STATUS "The installer will not relocate ${_openms_app_bundles_text} (${APPLICATIONS_COMPONENT_PLIST})")

unset(_openms_app_bundles)
unset(_openms_app_bundles_text)
unset(_openms_bundle_dir)
unset(_openms_plist_entries)

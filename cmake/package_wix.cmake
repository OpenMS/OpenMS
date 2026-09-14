# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer, Chris Bielow $
# --------------------------------------------------------------------------

## Windows installer -- WiX/MSI variant (prototype).
##
## This is an alternative to cmake/package_nsis.cmake. It produces a Windows
## Installer (.msi) database via CPack's WIX generator instead of an NSIS .exe.
##
## Scope of this prototype -- deliberately the *minimal* path:
##   * perMachine install only (NSIS today offers a per-user/per-machine choice
##     via the NsisMultiUser plugin; CPACK_WIX_INSTALL_SCOPE accepts perMachine
##     or perUser but has no combined mode, so that choice is not reproducible)
##   * English only (the NSIS installer offers German/English/French)
##   * no splash screen (no MSI equivalent; would need a Burn bootstrapper)
##   * no migration from an existing NSIS-based installation (see NOTE below)
##
## What comes for free relative to the NSIS template, and is therefore *not*
## reimplemented here:
##   * uninstall bookkeeping        (AdvUninstLog2.nsh -- MSI tracks files itself)
##   * previous-version removal     (<MajorUpgrade> is in CPack's default template)
##   * reboot-on-locked-file, empty directory cleanup
##   * PATH edits longer than 1k    (no need for the patched 8k build of NSIS)

if("${CPACK_GENERATOR}" STREQUAL "WIX")
  ## CPACK_WIX_VERSION 4 selects the .NET 'wix' tool, which covers WiX v4 and v5.
  ## WiX v3 (candle/light) is EOL and needs .NET Framework 3.5; don't use it.
  set(CPACK_WIX_VERSION 4)

  find_program(WIX_EXE wix)
  if(NOT WIX_EXE)
    message(FATAL_ERROR
      "Could not find the 'wix' executable. Install it with\n"
      "  dotnet tool install --global wix --version 5.0.2\n"
      "and make sure %USERPROFILE%\\.dotnet\\tools is in $PATH.")
  endif()

  ## The WixUI dialog set lives in an extension that must be installed
  ## separately; 'wix build' fails late and unhelpfully if it is missing.
  execute_process(COMMAND ${WIX_EXE} extension list --global
                  OUTPUT_VARIABLE WIX_EXTENSION_LIST
                  ERROR_QUIET)
  if(NOT WIX_EXTENSION_LIST MATCHES "WixToolset\\.UI\\.wixext")
    message(FATAL_ERROR
      "The WiX UI extension is not installed. Install it with\n"
      "  wix extension add --global WixToolset.UI.wixext/5.0.2")
  endif()
endif()

## VS2019+ requires the architecture via -A / CMAKE_GENERATOR_PLATFORM rather
## than a generator name suffix, same as package_nsis.cmake.
if (CMAKE_VS_PLATFORM_NAME MATCHES ".*Win32.*" OR CMAKE_GENERATOR MATCHES ".*Win32.*")
  set(PLATFORM "32")
  set(CPACK_WIX_ARCHITECTURE "x86")
else()
  set(PLATFORM "64")
  set(CPACK_WIX_ARCHITECTURE "x64")
endif()

#### Install MSVC runtime libraries into /bin so we do not have to ship a
#### VC-Redist.exe alongside the installer (same as the NSIS path).
set(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE)
set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION ${INSTALL_LIB_DIR})
include(InstallRequiredSystemLibraries)

## ---------------------------------------------------------------------------
## Identity
## ---------------------------------------------------------------------------

## The UpgradeCode identifies "OpenMS" across all versions and MUST stay fixed
## forever: it is what lets an MSI recognise and replace an older OpenMS.
## Changing it makes new installers sit side by side with old ones instead of
## upgrading them. Generated once for this prototype; do not regenerate.
set(CPACK_WIX_UPGRADE_GUID "6F2B4E5A-9C3D-4A17-8E52-0B7D1A4C6F38")

## MSI ProductVersion only compares major.minor.build and caps them at
## 255.255.65535; anything beyond that (including a pre-release suffix) is
## either rejected or silently ignored during upgrade detection.
if(NOT OPENMS_PACKAGE_VERSION MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)$")
  message(FATAL_ERROR
    "OPENMS_PACKAGE_VERSION ('${OPENMS_PACKAGE_VERSION}') is not a plain "
    "major.minor.patch triple; MSI cannot express it as a ProductVersion.")
endif()
if(CMAKE_MATCH_1 GREATER 255 OR CMAKE_MATCH_2 GREATER 255 OR CMAKE_MATCH_3 GREATER 65535)
  message(FATAL_ERROR
    "OPENMS_PACKAGE_VERSION ('${OPENMS_PACKAGE_VERSION}') exceeds the MSI "
    "ProductVersion limits of 255.255.65535.")
endif()

## NOTE: nightly/pre-release builds carry a suffix in
## OPENMS_PACKAGE_VERSION_FULLSTRING (e.g. '3.5.0-pre-...'). That suffix can
## live in the *file name* but not in the MSI ProductVersion, so two different
## nightlies of the same base version are indistinguishable to the installer
## and upgrade over each other. AllowSameVersionUpgrades in CPack's default
## template makes that work, but "which nightly is installed" is no longer
## answerable from Add/Remove Programs alone.

if((DEFINED ENV{CPACK_PACKAGE_FILE_NAME}) AND (NOT "$ENV{CPACK_PACKAGE_FILE_NAME}" STREQUAL ""))
  set(CPACK_PACKAGE_FILE_NAME "$ENV{CPACK_PACKAGE_FILE_NAME}")
else()
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION_FULLSTRING}-Win${PLATFORM}")
endif()

## ---------------------------------------------------------------------------
## Install scope / UI
## ---------------------------------------------------------------------------

set(CPACK_WIX_INSTALL_SCOPE "perMachine")

## Our own shortcuts are authored in openms_extras.wxs; '.' stops CPack from
## additionally generating its own start menu folder and uninstall shortcut.
set(CPACK_WIX_PROGRAM_MENU_FOLDER ".")

set(CPACK_WIX_PRODUCT_ICON "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS.ico")

## NOTE: no CPACK_WIX_UI_BANNER / CPACK_WIX_UI_DIALOG here. WiX requires exactly
## 493x58 and 493x312; the existing NSIS artwork in cmake/Windows/images is
## 150x57 and 164x314, so it cannot be reused and new artwork has to be drawn.

## License.txt is plain text; the WIX generator converts .txt to RTF itself.
## (CPACK_RESOURCE_FILE_LICENSE is set in package_general.cmake.)

## ---------------------------------------------------------------------------
## Monolithic vs. component install
## ---------------------------------------------------------------------------
## Kept monolithic for this prototype, exactly as package_nsis.cmake does.
## Setting OPENMS_WIX_COMPONENTS=ON exposes CPack's component graph as an MSI
## feature tree -- see cmake/package_components.cmake and the notes in
## cmake/Windows/WiX/README.md; the component names are not currently
## consistent enough for that to produce a sensible feature tree.
option(OPENMS_WIX_COMPONENTS "Expose CPack components as an MSI feature tree (experimental)" OFF)
if(NOT OPENMS_WIX_COMPONENTS)
  set(CPACK_COMPONENT_ALL_IN_ONE 1)
  set(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
  set(CPACK_MONOLITHIC_INSTALL 1)
  set(CPACK_WIX_UI_REF "WixUI_InstallDir")
else()
  set(CPACK_WIX_UI_REF "WixUI_FeatureTree")
endif()

## ---------------------------------------------------------------------------
## Hand-authored fragments (PATH, shortcuts, file associations)
## ---------------------------------------------------------------------------

set(OPENMS_WIX_REGKEY "Software\\OpenMS\\OpenMS")
set(OPENMS_WIX_MENU_FOLDER "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION}")
## MSI directory properties resolve with a trailing backslash, so
## "[INSTALL_ROOT]bin" is the bin directory. Our INSTALL_*_DIR variables use
## forward slashes; WiX/MSI wants backslashes inside Formatted strings.
string(REPLACE "/" "\\" OPENMS_WIX_BIN_DIR_W "${INSTALL_BIN_DIR}")
string(REPLACE "/" "\\" OPENMS_WIX_SHARE_DIR_W "${INSTALL_SHARE_DIR}")
## Environment/@System must track the install scope.
if(CPACK_WIX_INSTALL_SCOPE STREQUAL "perMachine")
  set(OPENMS_WIX_ENV_SYSTEM "yes")
else()
  set(OPENMS_WIX_ENV_SYSTEM "no")
endif()

## PATH entries for each installed third party tool directory. The NSIS script
## discovers these at *install* time by walking share/OpenMS/THIRDPARTY
## (NSIS.template.in:471). MSI's Environment table is a build-time table, so
## the same list has to be resolved now -- which is fine, because
## install_thirdparty_folder() already collected it in THIRDPARTY_COMPONENT_GROUP.
set(OPENMS_WIX_THIRDPARTY_PATH_COMPONENTS "")
foreach(_tp IN LISTS THIRDPARTY_COMPONENT_GROUP)
  ## Component/Environment ids must be valid MSI identifiers.
  string(MAKE_C_IDENTIFIER "${_tp}" _tp_id)
  string(APPEND OPENMS_WIX_THIRDPARTY_PATH_COMPONENTS
"      <Component Id=\"OpenMS_Path_TP_${_tp_id}\">
        <RegistryValue Root=\"HKMU\" Key=\"${OPENMS_WIX_REGKEY}\\Components\"
                       Name=\"PathTP_${_tp_id}\" Type=\"integer\" Value=\"1\" KeyPath=\"yes\"/>
        <Environment Id=\"OpenMSPathTP_${_tp_id}\" Name=\"PATH\"
                     Value=\"[INSTALL_ROOT]${OPENMS_WIX_SHARE_DIR_W}\\THIRDPARTY\\${_tp}\"
                     Action=\"set\" Part=\"last\" Permanent=\"no\" System=\"${OPENMS_WIX_ENV_SYSTEM}\"/>
      </Component>
")
endforeach()

## File type registration, mirroring the OpenMSGUIExtensions macro in
## NSIS.template.in:69.
## NOTE: these variable names are spelled with the exact case of the executable
## because they are looked up as OPENMS_WIX_ASSOC_${_app} below, and CMake
## variable names are case sensitive.
set(OPENMS_WIX_ASSOC_TOPPView mzData mzXML mzML sqMass dta dta2D cdf idXML featureXML consensusXML)
set(OPENMS_WIX_ASSOC_TOPPAS   toppas)

set(OPENMS_WIX_FILE_ASSOCIATION_COMPONENTS "")
foreach(_app TOPPView TOPPAS)
  if(NOT DEFINED OPENMS_WIX_ASSOC_${_app})
    message(FATAL_ERROR "No file extensions declared for ${_app} (OPENMS_WIX_ASSOC_${_app} is unset).")
  endif()
  foreach(_ext IN LISTS OPENMS_WIX_ASSOC_${_app})
    string(MAKE_C_IDENTIFIER "${_ext}" _ext_id)
    string(APPEND OPENMS_WIX_FILE_ASSOCIATION_COMPONENTS
"      <Component Id=\"OpenMS_Assoc_${_ext_id}\">
        <RegistryValue Root=\"HKMU\" Key=\"Software\\Classes\\.${_ext}\"
                       Value=\"OpenMS.${_ext}\" Type=\"string\" KeyPath=\"yes\"/>
        <RegistryValue Root=\"HKMU\" Key=\"Software\\Classes\\OpenMS.${_ext}\"
                       Value=\"OpenMS ${_ext} file\" Type=\"string\"/>
        <RegistryValue Root=\"HKMU\" Key=\"Software\\Classes\\OpenMS.${_ext}\\DefaultIcon\"
                       Value=\"[INSTALL_ROOT]${OPENMS_WIX_SHARE_DIR_W}\\OpenMS_${_app}.ico\" Type=\"string\"/>
        <RegistryValue Root=\"HKMU\" Key=\"Software\\Classes\\OpenMS.${_ext}\\shell\\open\\command\"
                       Value=\"&quot;[INSTALL_ROOT]${OPENMS_WIX_BIN_DIR_W}\\${_app}.exe&quot; &quot;%1&quot;\" Type=\"string\"/>
      </Component>
")
  endforeach()
endforeach()

configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/Windows/WiX/openms_extras.wxs.in"
  "${PROJECT_BINARY_DIR}/openms_extras.wxs"
  @ONLY)

set(CPACK_WIX_EXTRA_SOURCES "${PROJECT_BINARY_DIR}/openms_extras.wxs")
set(CPACK_WIX_PATCH_FILE   "${PROJECT_SOURCE_DIR}/cmake/Windows/WiX/openms_patch.xml")
set(CPACK_WIX_EXTENSIONS   "WixToolset.UI.wixext")

## The .ico files used by the file associations above are installed by the NSIS
## script directly out of the source tree; with WiX they have to be part of the
## normal install payload so that MSI tracks and removes them.
install(FILES
          "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS_TOPPView.ico"
          "${PROJECT_SOURCE_DIR}/cmake/Windows/OpenMS_TOPPAS.ico"
        DESTINATION ${INSTALL_SHARE_DIR}
        COMPONENT share)

## Likewise the release notes, which the NSIS "-License" section pulls straight
## from cmake/Windows.
install(FILES "${PROJECT_SOURCE_DIR}/cmake/Windows/ReleaseNotes.txt"
        DESTINATION .
        COMPONENT share)

## ---------------------------------------------------------------------------
## Targets
## ---------------------------------------------------------------------------

add_custom_target(dist
  COMMAND cpack -G ${CPACK_GENERATOR}
  COMMENT "Building ${CPACK_GENERATOR} package"
)

if (DEFINED SIGNING_IDENTITY AND NOT "${SIGNING_IDENTITY}" STREQUAL "")
  add_custom_target(signed_dist
    COMMAND signtool sign /v /n "${SIGNING_IDENTITY}" /t http://timestamp.digicert.com ${CPACK_PACKAGE_FILE_NAME}.msi
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Signing ${CPACK_PACKAGE_FILE_NAME}.msi with '${SIGNING_IDENTITY}'"
    DEPENDS dist
  )
endif()

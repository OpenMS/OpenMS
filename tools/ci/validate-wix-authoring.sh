#!/usr/bin/env bash
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
#
# Lint the hand-authored WiX sources (cmake/Windows/WiX/) on any platform.
#
# Building an actual .msi requires Windows: the WiX toolset says so itself on
# startup, and the final bind step P/Invokes msi.dll to create the installer
# database. Everything *before* that step -- XML parsing, schema validation,
# reference resolution and linking -- is pure managed code and runs fine on
# Linux. 'wix build -o <file>.wixipl' stops exactly at that boundary, which
# makes it a usable authoring linter off Windows.
#
# What this catches: malformed XML, misspelled elements/attributes, component
# groups that the patch file references but that do not exist (which would
# otherwise be silently dropped from the .msi), and broken CMake substitution
# into openms_extras.wxs.in.
#
# What it does NOT catch: anything ICE validation, installer UI, or actually
# running the installer would find. Those still need a Windows runner.
#
# Usage: tools/ci/validate-wix-authoring.sh

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
work_dir="$(mktemp -d)"
trap 'rm -rf "${work_dir}"' EXIT

WIX_VERSION="5.0.2"
dotnet_root="${DOTNET_ROOT:-/usr/share/dotnet}"

## Reuse the repo's pinned, checksum-verified SDK installer rather than
## dotnet-install.sh (see the rationale in that script's header).
if [[ ! -x "${dotnet_root}/dotnet" ]]; then
  "${repo_root}/tools/ci/install_dotnet_sdk_linux.sh" "${dotnet_root}"
fi
export DOTNET_ROOT="${dotnet_root}"
export DOTNET_CLI_TELEMETRY_OPTOUT=1
export PATH="${dotnet_root}:${HOME}/.dotnet/tools:${PATH}"

if ! command -v wix >/dev/null 2>&1; then
  dotnet tool install --global wix --version "${WIX_VERSION}" >/dev/null
fi
wix extension add --global "WixToolset.UI.wixext/${WIX_VERSION}" >/dev/null 2>&1 || true

## ---------------------------------------------------------------------------
## 1. Generate openms_extras.wxs the same way a real configure would.
## ---------------------------------------------------------------------------
mkdir -p "${work_dir}/fixture"
cat > "${work_dir}/fixture/CMakeLists.txt" <<'FIXEOF'
cmake_minimum_required(VERSION 3.24)
project(openms_wix_lint NONE)

# Stand-ins for the OpenMS variables that cmake/package_wix.cmake reads.
set(PROJECT_SOURCE_DIR "${OPENMS_ROOT}")
set(INSTALL_BIN_DIR bin)
set(INSTALL_LIB_DIR bin)
set(INSTALL_SHARE_DIR share/OpenMS)
set(OPENMS_PACKAGE_VERSION "3.5.0")
set(OPENMS_PACKAGE_VERSION_FULLSTRING "3.5.0")
set(CPACK_PACKAGE_NAME "OpenMS")
# A representative sample of what install_thirdparty_folder() collects.
set(THIRDPARTY_COMPONENT_GROUP Comet Crux MSGFPlus Sirius XTandem pwiz-bin)
# Left empty so package_wix.cmake skips its wix/extension discovery block;
# this lint does not need the real toolchain lookup.
set(CPACK_GENERATOR "")

include("${OPENMS_ROOT}/cmake/package_wix.cmake")
FIXEOF

cmake -S "${work_dir}/fixture" -B "${work_dir}/fixture-build" \
      -D "OPENMS_ROOT=${repo_root}" > "${work_dir}/configure.log" 2>&1 || {
  echo "FAIL: configuring the lint fixture failed" >&2
  cat "${work_dir}/configure.log" >&2
  exit 1
}

extras="${work_dir}/fixture-build/openms_extras.wxs"
if [[ ! -f "${extras}" ]]; then
  echo "FAIL: package_wix.cmake did not produce openms_extras.wxs" >&2
  exit 1
fi

## configure_file() substitutes an *undefined* @VAR@ with the empty string
## rather than leaving the placeholder in place, so a mistyped variable name in
## the .wxs.in does not show up as leftover @...@ -- it shows up as an empty
## attribute value (e.g. Name="") that WiX is happy to accept. Check for both:
## the placeholder form catches a literal that was never meant to be a
## placeholder, the empty-attribute form catches the typo.
if grep -n '@[A-Z_]\{4,\}@' "${extras}"; then
  echo "FAIL: unsubstituted placeholders remain in openms_extras.wxs" >&2
  exit 1
fi
if grep -n '=""' "${extras}"; then
  echo "FAIL: empty attribute value in openms_extras.wxs" >&2
  echo "      (an undefined variable in openms_extras.wxs.in substitutes to nothing)" >&2
  exit 1
fi

## Content checks. Linking alone does not catch a fragment that came out empty:
## an empty <ComponentGroup> is perfectly valid WiX and links without complaint,
## so a mis-spelled CMake variable feeding one of the foreach() loops would
## silently ship an installer with, say, no file associations at all. Assert the
## expected counts instead.
check_count() {
  local label="$1" pattern="$2" expected="$3"
  local actual
  actual="$(grep -c "${pattern}" "${extras}" || true)"
  if [[ "${actual}" -ne "${expected}" ]]; then
    echo "FAIL: expected ${expected} ${label} in openms_extras.wxs, found ${actual}" >&2
    return 1
  fi
  echo "  ${label}: ${actual}"
}
## 1 for bin/ plus one per third party tool in the fixture's
## THIRDPARTY_COMPONENT_GROUP (6 entries).
check_count "PATH components"        'Id="OpenMS_Path_'  7
## 10 TOPPView extensions + 1 TOPPAS extension, see OPENMS_WIX_ASSOC_* .
check_count "file associations"      'Id="OpenMS_Assoc_' 11
check_count "start menu shortcuts"   '<Shortcut '        4

## ---------------------------------------------------------------------------
## 2. Stand in for the sources CPack's WIX generator emits at packaging time.
## ---------------------------------------------------------------------------
## The real build gets INSTALL_ROOT and <Feature Id="ProductFeature"> from
## CPack, with our ComponentGroupRefs injected via CPACK_WIX_PATCH_FILE. Here
## they are written out directly so the fragments have something to link
## against. The ComponentGroupRef list is read out of the real patch file, so
## a group renamed in only one of the two files is caught.
refs="$(grep -o 'ComponentGroupRef Id="[^"]*"' \
        "${repo_root}/cmake/Windows/WiX/openms_patch.xml" \
        | sed 's/.*Id="\([^"]*\)"/      <ComponentGroupRef Id="\1"\/>/')"
if [[ -z "${refs}" ]]; then
  echo "FAIL: no ComponentGroupRef entries found in openms_patch.xml" >&2
  exit 1
fi

cat > "${work_dir}/stub_package.wxs" <<STUBEOF
<?xml version="1.0" encoding="UTF-8"?>
<Wix xmlns="http://wixtoolset.org/schemas/v4/wxs">
  <Package Name="OpenMS" Version="3.5.0" Manufacturer="OpenMS.de"
           UpgradeCode="6F2B4E5A-9C3D-4A17-8E52-0B7D1A4C6F38"
           Scope="perMachine" Compressed="yes">
    <MajorUpgrade DowngradeErrorMessage="A later version is already installed."/>
    <StandardDirectory Id="ProgramFiles64Folder">
      <Directory Id="INSTALL_ROOT"/>
    </StandardDirectory>
    <Feature Id="ProductFeature">
${refs}
    </Feature>
  </Package>
</Wix>
STUBEOF

## ---------------------------------------------------------------------------
## 3. Lint.
## ---------------------------------------------------------------------------
## Directory/@Name is dropped first: on Linux *every* value of that attribute
## trips "WIX0389: ... is not a relative path", because the check runs through
## Windows path semantics. The attribute only supplies the on-disk folder name
## and has no bearing on the reference resolution this lint is for, so removing
## it keeps the rest of the file checkable. Verified against wix 5.0.2.
sed 's/\(<Directory[^>]*\) Name="[^"]*"/\1/' "${extras}" > "${work_dir}/extras_lint.wxs"

set +e
wix build "${work_dir}/extras_lint.wxs" "${work_dir}/stub_package.wxs" \
          -ext WixToolset.UI.wixext \
          -o "${work_dir}/out.wixipl" > "${work_dir}/wix.log" 2>&1
status=$?
set -e

grep -v "only supports Windows" "${work_dir}/wix.log" | grep -v '^[[:space:]]*$' || true

if [[ ${status} -ne 0 || ! -f "${work_dir}/out.wixipl" ]]; then
  echo "FAIL: WiX authoring did not compile and link" >&2
  exit 1
fi

echo "OK: WiX authoring compiles and links (msi generation still needs Windows)"

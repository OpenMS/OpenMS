#!/usr/bin/env bash
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
#
# Fetch the assets the pyOpenMS wheel builds need for native Thermo RAW support,
# without going through NuGet or the .NET SDK's package restore:
#
#   1. the pre-built managed half of openms-thermo-bridge (ThermoWrapperManaged.dll,
#      its runtimeconfig.json and the Thermo CommonCore assemblies), published as a
#      zip on the bridge's GitHub release. Passing the extracted directory to CMake
#      as OPENMS_THERMO_BRIDGE_PREBUILT_MANAGED_DIR skips the bridge's own
#      'dotnet publish' step, and with it the five Thermo .nupkg downloads from
#      raw.githubusercontent.com and the nuget.org restore of their dependencies.
#      Compiling the native bridge still needs the nethost headers from a .NET
#      SDK / host pack.
#   2. PNNL's Angiotensin_AllScans.raw from archive.openms.de, used to smoke-test
#      the finished wheel.
#
# Every download is retried and verified against a pinned SHA-256, and the bridge
# tag is checked against the FetchContent pin in cmake/cmake_findExternalLibs.cmake
# so the two cannot drift apart silently.
#
# Usage: tools/ci/fetch_thermo_assets.sh <linux-x64|osx-arm64|win-x64> [dest-dir]
#   dest-dir defaults to the current directory; creates <dest-dir>/thermo-managed
#   and <dest-dir>/thermo-testdata. Existing, verified files are reused, so the
#   two directories can be put into a CI cache.

set -euo pipefail

BRIDGE_TAG="v0.2.3"
RAW_URL="https://archive.openms.de/openms/testfiles/Angiotensin_AllScans.raw"
RAW_SHA256="3a0236f719e7c91e3c958f57f4e66ae422803ec3e6a997b9af4d2af395332b9f"

platform="${1:-}"
dest="${2:-.}"
# Plain case instead of an associative array: the stock macOS Bash is 3.2, which
# has no 'declare -A', and the test-wheels job runs this script with that Bash.
case "$platform" in
  linux-x64) managed_sha256="c3388188f350280e1532d69e460c2cd76f5c69f4b4e6e4616a1a21e4347f5e4f" ;;
  osx-arm64) managed_sha256="0253e156630db5bca4cdd4d401cfec5bfc532349b73b48ed2a60a004e3ff4f80" ;;
  win-x64)   managed_sha256="ac737effe6e5c4b379bf9fe242ac7c969f7b61616741bfb7241f3fdc9c620581" ;;
  *)
    echo "usage: $0 <linux-x64|osx-arm64|win-x64> [dest-dir]" >&2
    exit 2
    ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"

# Guard against the pin in this script and the FetchContent pin diverging.
if ! grep -A6 'OpenMSThermoBridge$' "${repo_root}/cmake/cmake_findExternalLibs.cmake" \
     | grep -qE "GIT_TAG[[:space:]]+${BRIDGE_TAG}([[:space:]]|$)"; then
  echo "error: BRIDGE_TAG=${BRIDGE_TAG} in $0 does not match the GIT_TAG pinned for" >&2
  echo "       OpenMSThermoBridge in cmake/cmake_findExternalLibs.cmake" >&2
  exit 1
fi

# Hash via stdin: given a file *name*, GNU coreutils escapes backslashes in it and
# prefixes the whole line with '\' (seen with Git Bash on Windows, where
# $GITHUB_WORKSPACE is D:\a\...), which would corrupt the extracted digest.
sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum < "$1" | cut -d' ' -f1
  else
    shasum -a 256 < "$1" | cut -d' ' -f1
  fi
}

# download <url> <file> <sha256>: reuse a verified existing file, otherwise fetch
# with retries and fail loudly on a hash mismatch (never fall back silently).
download() {
  local url="$1" file="$2" expected="$3"
  if [[ -f "$file" ]] && [[ "$(sha256_of "$file")" == "$expected" ]]; then
    echo "reusing verified $file"
    return 0
  fi
  echo "downloading $url"
  rm -f "$file"
  curl --fail --location --silent --show-error \
       --retry 8 --retry-delay 5 --retry-all-errors --connect-timeout 30 --max-time 600 \
       --output "$file" "$url"
  local actual
  actual="$(sha256_of "$file")"
  if [[ "$actual" != "$expected" ]]; then
    echo "error: SHA-256 mismatch for $file" >&2
    echo "       expected $expected" >&2
    echo "       actual   $actual" >&2
    rm -f "$file"
    exit 1
  fi
}

mkdir -p "${dest}/thermo-managed" "${dest}/thermo-testdata"

# The managed directory is handed to CMake verbatim and installed as-is, so keep the
# zip out of it. A fully extracted directory (e.g. restored from a CI cache) is reused.
required_managed=(ThermoWrapperManaged.dll ThermoWrapperManaged.runtimeconfig.json ThermoFisher.CommonCore.RawFileReader.dll)
managed_complete=true
for required in "${required_managed[@]}"; do
  [[ -f "${dest}/thermo-managed/${required}" ]] || managed_complete=false
done
if [[ "$managed_complete" == true ]]; then
  echo "reusing extracted managed bridge in ${dest}/thermo-managed"
else
  zip_name="openms-thermo-bridge-managed-${platform}-${BRIDGE_TAG}.zip"
  zip_path="$(mktemp -d)/${zip_name}"
  download "https://github.com/jpfeuffer/openms-thermo-bridge/releases/download/${BRIDGE_TAG}/${zip_name}" \
           "$zip_path" "${managed_sha256}"
  # cmake -E tar is available on every runner and understands zip on all platforms.
  (cd "${dest}/thermo-managed" && cmake -E tar xf "$zip_path")
  rm -rf "$(dirname "$zip_path")"
  for required in "${required_managed[@]}"; do
    if [[ ! -f "${dest}/thermo-managed/${required}" ]]; then
      echo "error: ${required} missing after extracting ${zip_name}" >&2
      exit 1
    fi
  done
fi

download "$RAW_URL" "${dest}/thermo-testdata/Angiotensin_AllScans.raw" "$RAW_SHA256"

echo "Thermo assets ready:"
echo "  managed bridge: ${dest}/thermo-managed"
echo "  RAW test file:  ${dest}/thermo-testdata/Angiotensin_AllScans.raw"

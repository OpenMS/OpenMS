#!/usr/bin/env bash
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
#
# Install a pinned .NET SDK (linux-x64) from Microsoft's release archive into a
# directory, verifying the tarball against the SHA-512 Microsoft publishes in
# https://builds.dotnet.microsoft.com/dotnet/release-metadata/8.0/releases.json.
#
# The pyOpenMS wheel build needs this for exactly one thing: the nethost headers
# and static library in packs/Microsoft.NETCore.App.Host.linux-x64, which the
# native openms-thermo-bridge compiles against. It is used instead of
#   * the distro dotnet-sdk RPM: AlmaLinux names its host pack
#     Microsoft.NETCore.App.Host.rhel.9-x64, which the bridge's
#     FindDotNetHost.cmake (linux-x64 RID) does not find;
#   * dotnet-install.sh: a mutable script fetched from a redirecting endpoint
#     and executed unverified.
#
# Usage: tools/ci/install_dotnet_sdk_linux.sh [install-dir]   (default /usr/share/dotnet)
#   /usr/share/dotnet is a default search location for both FindDotNetHost.cmake
#   and the nethost runtime resolver, so DOTNET_ROOT is optional there.

set -euo pipefail

DOTNET_SDK_VERSION="8.0.424"
DOTNET_SDK_URL="https://builds.dotnet.microsoft.com/dotnet/Sdk/${DOTNET_SDK_VERSION}/dotnet-sdk-${DOTNET_SDK_VERSION}-linux-x64.tar.gz"
DOTNET_SDK_SHA512="6503fd9f464d5e3a4f43a881d2b74afc6a2c46ceda74d027f1565b7239f4b3ec884857c03c0dcd49eb52f384d5ae1fa5aaf135f0a6aabc5518103aceed643c74"

install_dir="${1:-/usr/share/dotnet}"

if [[ -f "${install_dir}/sdk/${DOTNET_SDK_VERSION}/dotnet.dll" ]]; then
  echo "reusing .NET SDK ${DOTNET_SDK_VERSION} in ${install_dir}"
else
  tmp="$(mktemp -d)"
  trap 'rm -rf "$tmp"' EXIT
  echo "downloading ${DOTNET_SDK_URL}"
  curl --fail --location --proto '=https' --proto-redir '=https' --silent --show-error \
       --retry 8 --retry-delay 5 --retry-all-errors --connect-timeout 30 --max-time 900 \
       --output "${tmp}/dotnet-sdk.tar.gz" "${DOTNET_SDK_URL}"
  actual="$(sha512sum < "${tmp}/dotnet-sdk.tar.gz" | cut -d' ' -f1)"
  if [[ "${actual}" != "${DOTNET_SDK_SHA512}" ]]; then
    echo "error: SHA-512 mismatch for dotnet-sdk-${DOTNET_SDK_VERSION}-linux-x64.tar.gz" >&2
    echo "       expected ${DOTNET_SDK_SHA512}" >&2
    echo "       actual   ${actual}" >&2
    exit 1
  fi
  mkdir -p "${install_dir}"
  tar -xzf "${tmp}/dotnet-sdk.tar.gz" -C "${install_dir}"
fi

host_pack="$(ls -d "${install_dir}"/packs/Microsoft.NETCore.App.Host.linux-x64/*/runtimes/linux-x64/native 2>/dev/null | head -1)"
if [[ -z "${host_pack}" || ! -f "${host_pack}/nethost.h" ]]; then
  echo "error: nethost host pack not found under ${install_dir}/packs" >&2
  exit 1
fi
echo ".NET SDK ${DOTNET_SDK_VERSION} ready in ${install_dir}; nethost pack: ${host_pack}"

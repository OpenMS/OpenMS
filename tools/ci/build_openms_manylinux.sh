#!/usr/bin/env bash
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
#
# Build and install the OpenMS core library (libOpenMS, libOpenSwathAlgo) that the
# Linux pyOpenMS wheels are built against. Runs as root inside the stock
# quay.io/pypa/manylinux_2_34_{x86_64,aarch64} image, with the repository mounted at
# the directory given as first argument (default /project); the wheel workflow calls
# it from `docker run` before cibuildwheel runs in the same image.
#
# The dependencies come from vcpkg (the repository's vcpkg submodule and manifest),
# built as static libraries with the overlay triplets x64-linux / arm64-linux, so
# the wheel links them into libOpenMS.so and auditwheel grafts only OpenMS's own
# shared libraries plus system libraries such as libgomp. Set VCPKG_BINARY_SOURCES
# to reuse binary packages of earlier runs.
#
# Results, all below the repository directory so cibuildwheel's copy of the project
# carries them:
#   install/                          the OpenMS installation (OpenMS_ROOT)
#   vcpkg_installed/<triplet>/        the vcpkg dependencies (CMAKE_PREFIX_PATH)
#
# Usage: tools/ci/build_openms_manylinux.sh [repository-dir]

set -euo pipefail

src="${1:-/project}"
cd "${src}"

case "$(uname -m)" in
  x86_64|amd64) triplet="x64-linux" ;;
  aarch64|arm64) triplet="arm64-linux" ;;
  *)
    echo "error: no vcpkg triplet for machine $(uname -m)" >&2
    exit 2
    ;;
esac

# The bind-mounted workspace is owned by another user, so without this git exits
# 128 and the wheels carry a placeholder version instead of the commit.
git config --global --add safe.directory "${src}"
git config --global --add safe.directory "${src}/vcpkg"

# Tools the vcpkg ports need on top of what the manylinux image has (gcc-toolset
# with gfortran for lapack-reference, CMake, git, curl, unzip, tar, pkg-config,
# bison and recent autotools in /usr/local/bin):
#   ninja-build, zip           vcpkg itself
#   autoconf-archive           the coin-or ports (autoreconf needs its AX_ macros)
#   flex                       thrift (for Arrow's Parquet support)
#   perl modules               OpenSSL's Configure
#   libicu                     the .NET SDK installed below
dnf install -y \
  ninja-build \
  zip \
  autoconf-archive \
  flex \
  perl-IPC-Cmd \
  perl-FindBin \
  perl-File-Compare \
  perl-File-Copy \
  perl-Time-Piece \
  perl-lib \
  libicu

# The openms-thermo-bridge port compiles against the nethost headers of a .NET SDK
# and publishes the managed assemblies with it. The pinned, SHA-512 verified SDK
# tarball is used instead of the distro dotnet-sdk-8.0 RPM, whose host pack is
# named with the rhel.9 RID that the bridge does not search (see the script).
bash tools/ci/install_dotnet_sdk_linux.sh /usr/share/dotnet
export DOTNET_ROOT=/usr/share/dotnet
export PATH="${DOTNET_ROOT}:${PATH}"
export DOTNET_CLI_TELEMETRY_OPTOUT=1
export DOTNET_NOLOGO=1

./vcpkg/bootstrap-vcpkg.sh -disableMetrics

# Manifest features: the options this build turns on that have a vcpkg port
# (WITH_THERMO_RAW and WITH_OPENTIMS, both ON by default). The libraries OpenMS
# vendors (USE_EXTERNAL_* OFF) stay vendored, as in the earlier contrib-based wheels.
#
# ARROW_USE_STATIC is ON by default; it is spelled out because the standalone
# pyOpenMS build that links _arrow_zerocopy against the same Arrow is configured
# with it too (CIBW_ENVIRONMENT in the workflow).
#
# --clean-after-build drops each port's build tree once it is installed, which
# keeps the cold build (Arrow with the AWS SDK, Boost, ...) within the runner's disk.
cmake -S . -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE="${src}/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DOPENMS_USE_VCPKG=ON \
  -DVCPKG_TARGET_TRIPLET="${triplet}" \
  -DVCPKG_HOST_TRIPLET="${triplet}" \
  -DVCPKG_INSTALLED_DIR="${src}/vcpkg_installed" \
  -DVCPKG_MANIFEST_FEATURES="openms-thermo-bridge;opentims" \
  -DVCPKG_INSTALL_OPTIONS="--clean-after-build" \
  -DARROW_USE_STATIC=ON \
  -DCMAKE_INSTALL_PREFIX="${src}/install" \
  -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
  -DPYOPENMS=OFF \
  -DWITH_GUI=OFF \
  -DBUILD_TOPP_TOOLS=OFF \
  -DINSTALL_OPENMS_EXAMPLES=OFF \
  -DWITH_THERMO_RAW=ON \
  -DENABLE_TUTORIALS=OFF \
  -DENABLE_DOCS=OFF

# core layer only (libOpenMS, libOpenSwathAlgo): the CLI layer is a separate install component
cmake --build build --target OpenMS OpenSwathAlgo
cmake --install build --component library --strip
cmake --install build --component OpenMS_headers
cmake --install build --component OpenSwathAlgo_headers
cmake --install build --component thirdparty_headers
cmake --install build --component share
cmake --install build --component cmake

# cibuildwheel copies the whole project directory into its container. The OpenMS
# build tree and vcpkg's download and package staging areas are not needed there.
rm -rf build vcpkg/buildtrees vcpkg/packages vcpkg/downloads

echo "OpenMS installed to ${src}/install; vcpkg dependencies in ${src}/vcpkg_installed/${triplet}"

#!/usr/bin/env bash

set -eu
set -o pipefail

# Code between the following doxygen markers are included in the
# public-facing OpenMS installation instructions.

# [installation_documentation]
# Add "universe" and update:
sudo add-apt-repository universe
sudo apt update

# CMake 4.1+ (required for C++23 import std with GCC)
# Remove any pre-installed CMake (GitHub runners ship 3.31.x in /usr/local/bin
# which shadows the Kitware version in /usr/bin)
sudo rm -f /usr/local/bin/cmake /usr/local/bin/ctest /usr/local/bin/cpack
wget -qO- https://apt.kitware.com/keys/kitware-archive-latest.asc | sudo tee /etc/apt/trusted.gpg.d/kitware.asc
echo "deb https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" | sudo tee /etc/apt/sources.list.d/kitware.list
sudo apt-get update

# Required dependencies (install before GCC 15 so noble packages resolve
# against the default noble runtime libraries):
sudo apt-get -qq install -y \
  build-essential \
  cmake \
  autoconf \
  patch \
  libtool \
  git \
  automake \
  ninja-build \
  xvfb \
  ccache \
  qt6-base-dev \
  libqt6svg6-dev \
  libqt6opengl6-dev \
  libqt6openglwidgets6 \
  libgl-dev \
  libeigen3-dev \
  libboost-random-dev \
  libboost-regex-dev \
  libboost-iostreams-dev \
  libboost-date-time-dev \
  libboost-math-dev \
  libxerces-c-dev \
  zlib1g-dev \
  libsvm-dev \
  libbz2-dev \
  coinor-libcoinmp-dev \
  libhdf5-dev \
  libsqlite3-dev \
  libsqlitecpp-dev \
  nlohmann-json3-dev \
  libsimde-dev

# GCC 15 for C++23 import std support
# Ubuntu 24.04 (noble) doesn't ship GCC 15. Install from Ubuntu 25.10
# (questing) which has GCC 15.2 with libstdc++.modules.json for import std.
# Plucky (25.04) only has GCC 15.0.1 (pre-release) which lacks module support.
# APT pinning ensures only GCC packages are pulled from questing.
# NOTE: GCC 15 is installed AFTER the other dependencies above because it
# upgrades shared runtime libraries (libstdc++6, libgcc-s1) from 14 to 15,
# which would cause dependency resolution failures for noble packages.
if [ "$(uname -m)" = "x86_64" ]; then
    GCC_REPO_URL="http://archive.ubuntu.com/ubuntu"
else
    GCC_REPO_URL="http://ports.ubuntu.com/ubuntu-ports"
fi
echo "deb $GCC_REPO_URL questing main universe" | sudo tee /etc/apt/sources.list.d/questing-gcc.list
# Pin questing packages low by default; only GCC 15 packages get elevated priority
sudo tee /etc/apt/preferences.d/questing-pin > /dev/null <<'PINEOF'
Package: *
Pin: release n=questing
Pin-Priority: 100

Package: gcc-15* g++-15* cpp-15* libgcc-15-dev libstdc++-15-dev libstdc++6 libgcc-s1 libcc1-0 libgomp1 libitm1 libatomic1 libasan8 liblsan0 libtsan2 libubsan1 libhwasan0 libquadmath0 binutils* libbinutils* libctf* libgprofng* libsframe*
Pin: release n=questing
Pin-Priority: 500
PINEOF
sudo apt-get update
sudo apt-get install -y g++-15 binutils
# Verify the compiler and assembler are accessible
g++-15 --version
as --version

# Install uv (Python package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install Apache Arrow
sudo apt-get install -y -V ca-certificates lsb-release wget
wget https://packages.apache.org/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
sudo apt-get install -y -V ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
sudo apt update
# Install libcurl-dev as a workaround for Arrow CMake config issue
# See: https://github.com/apache/arrow/issues/48885
sudo apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev \
      libarrow-dev \
      libparquet-dev

# Optional dependencies:
sudo apt-get -qq install -y \
  doxygen \
  ghostscript \
  graphviz
# [installation_documentation]

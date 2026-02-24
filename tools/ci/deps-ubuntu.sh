#!/usr/bin/env bash

set -eu
set -o pipefail

# Code between the following doxygen markers are included in the
# public-facing OpenMS installation instructions.

# [installation_documentation]
# Add "universe" and update:
sudo add-apt-repository universe
sudo apt update

# Required dependencies:
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

# Build minizip-ng from source (no apt package available on Ubuntu 24.04)
MINIZIP_NG_VERSION=4.0.7
wget -q "https://github.com/zlib-ng/minizip-ng/archive/refs/tags/${MINIZIP_NG_VERSION}.tar.gz" -O /tmp/minizip-ng-${MINIZIP_NG_VERSION}.tar.gz
tar xzf /tmp/minizip-ng-${MINIZIP_NG_VERSION}.tar.gz -C /tmp
cmake -S /tmp/minizip-ng-${MINIZIP_NG_VERSION} -B /tmp/minizip-ng-build \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON \
  -DMZ_FETCH_LIBS=OFF -DMZ_LIB_SUFFIX=-ng -DMZ_COMPAT=OFF \
  -DMZ_ZLIB=ON -DMZ_BZIP2=ON -DMZ_LZMA=OFF -DMZ_ZSTD=OFF \
  -DMZ_OPENSSL=OFF -DMZ_ICONV=OFF -DMZ_PKCRYPT=ON -DMZ_WZAES=OFF
cmake --build /tmp/minizip-ng-build -j$(nproc)
sudo cmake --install /tmp/minizip-ng-build

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

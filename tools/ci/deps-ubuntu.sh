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
  libhdf5-dev

  sudo apt update
  sudo apt-get install -y -V ca-certificates lsb-release wget
  wget https://packages.apache.org/artifactory/arrow/$(lsb_release --id --short | tr 'A-Z' 'a-z')/apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
  sudo apt update
  sudo apt-get install -y -V ./apache-arrow-apt-source-latest-$(lsb_release --codename --short).deb
  sudo apt update
  sudo apt-get install -y --no-install-recommends \
        libarrow-dev \
        libparquet-dev \

# Optional dependencies:
sudo apt-get -qq install -y \
  doxygen \
  ghostscript \
  graphviz

# Install HIGHS LP solver from source
# HiGHS is a modern LP/MIP solver that exports CMake config files
HIGHS_VERSION=v1.7.2
cd /tmp
git clone --depth 1 --branch ${HIGHS_VERSION} https://github.com/ERGO-Code/HiGHS.git
cd HiGHS
cmake -S. -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --build build --parallel
sudo cmake --install build
cd /
rm -rf /tmp/HiGHS
# [installation_documentation]

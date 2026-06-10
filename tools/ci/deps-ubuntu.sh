#!/usr/bin/env bash

set -eu
set -o pipefail

# Parse command line arguments
SKIP_GUI_DEPS=false

while [[ $# -gt 0 ]]; do
  case $1 in
    --skip-gui-deps)
      SKIP_GUI_DEPS=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 [--skip-gui-deps]"
      exit 1
      ;;
  esac
done

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

# GUI dependencies (can be skipped for non-GUI builds):
if [ "$SKIP_GUI_DEPS" = false ]; then
  sudo apt-get -qq install -y \
    qt6-base-dev \
    libqt6svg6-dev \
    libqt6opengl6-dev \
    libqt6openglwidgets6
fi

# libzip (ZIP64 archive support) — available via apt on Ubuntu 24.04
sudo apt-get -qq install -y libzip-dev

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

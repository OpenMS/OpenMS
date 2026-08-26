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
  autoconf \
  autoconf-archive \
  automake \
  build-essential \
  cmake \
  git \
  libtool \
  ninja-build \
  patch \
  pkg-config

# GUI dependencies (can be skipped for non-GUI builds):
if [ "$SKIP_GUI_DEPS" = false ]; then
  sudo apt-get -qq install -y \
    qt6-base-dev \
    libqt6svg6-dev \
    libqt6opengl6-dev \
    libqt6openglwidgets6
fi

# Optional dependencies:
sudo apt-get -qq install -y \
  doxygen \
  graphviz

# [installation_documentation]

# These are only needed in CI:
sudo apt-get -qq install -y \
  ccache \
  rclone

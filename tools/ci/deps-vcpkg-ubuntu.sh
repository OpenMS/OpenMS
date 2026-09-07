#!/usr/bin/env bash
set -e
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

# [installation_documentation]
sudo add-apt-repository universe
sudo apt-get update
sudo apt-get -qq install -y \
  build-essential \
  cmake \
  ninja-build \
  autoconf \
  autoconf-archive \
  automake \
  libtool \
  patch \
  ccache \
  curl \
  git \
  zip \
  unzip \
  pkg-config \
  flex \
  bison


# GUI dependencies (can be skipped for non-GUI builds):
if [ "${SKIP_GUI_DEPS:-false}" = false ]; then
  # GUI dependencies (libgl-dev, xvfb) for building/testing TOPPView, TOPPAS:
  sudo apt-get -qq install -y \
    libgl-dev \
    xvfb
    
  sudo apt-get -qq install -y \
    qt6-base-dev \
    libqt6svg6-dev \
    libqt6opengl6-dev \
    libqt6openglwidgets6
fi

# Install uv (Python package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Optional dependencies:
sudo apt-get -qq install -y \
  doxygen \
  graphviz
# [installation_documentation]
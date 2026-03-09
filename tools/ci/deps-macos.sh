#!/usr/bin/env bash

set -eu
set -o pipefail

# Unfortunately GitHub's macOS runner already has Python installed so
# we need to tell brew to overwrite the existing links.  The following
# function will be called when the brew commands below are executed.
# It then calls the real brew command.
function brew() {
  local action=$1
  shift

  # Bash on macOS doesn't allow using empty arrays.  Therefore we put
  # the action name in the flags array so it always has at least one
  # element.  This is also why we install bash below.
  local -a flags=("$action")

  if [ "$action" = "install" ]; then
    flags+=("--overwrite")
  fi

  command brew "${flags[@]}" "$@"
}

# Code between the following doxygen markers are included in the
# public-facing OpenMS installation instructions.

# [installation_documentation]
# Update the package index:
brew update

# Required dependencies:
brew install \
  python \
  ccache \
  autoconf \
  automake \
  libtool \
  ninja \
  libomp \
  libsvm \
  xerces-c \
  boost \
  sqlite \
  coinutils \
  cbc \
  cgl \
  clp \
  qtbase \
  qtsvg \
  apache-arrow \
  bash

# Install Eigen 3.4.0 from source (Homebrew's eigen is now 5.x which is incompatible)
# Use a separate install prefix to avoid conflicts with Homebrew
EIGEN_VERSION="3.4.0"
EIGEN_INSTALL_PREFIX="/usr/local"
EIGEN_SRC="/tmp/eigen-${EIGEN_VERSION}"

echo "Installing Eigen ${EIGEN_VERSION} from source..."

# Download Eigen
EIGEN_URL="https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz"
echo "Downloading from ${EIGEN_URL}..."
if ! curl -fSL --retry 3 --retry-delay 5 "${EIGEN_URL}" -o /tmp/eigen.tar.gz; then
  echo "Failed to download Eigen from GitLab, trying GitHub mirror..."
  curl -fSL --retry 3 "https://github.com/eigenteam/eigen-git-mirror/archive/refs/tags/${EIGEN_VERSION}.tar.gz" -o /tmp/eigen.tar.gz
  EIGEN_SRC="/tmp/eigen-git-mirror-${EIGEN_VERSION}"
fi

echo "Extracting Eigen..."
tar -xzf /tmp/eigen.tar.gz -C /tmp

echo "Listing /tmp for debug:"
ls -la /tmp/ | grep -i eigen || echo "No eigen directories found!"

if [ ! -d "${EIGEN_SRC}" ]; then
  echo "Expected directory ${EIGEN_SRC} not found, searching..."
  EIGEN_SRC=$(find /tmp -maxdepth 1 -type d -name "*eigen*" 2>/dev/null | head -1)
fi
echo "Using Eigen source at: ${EIGEN_SRC}"

if [ -z "${EIGEN_SRC}" ] || [ ! -d "${EIGEN_SRC}" ]; then
  echo "ERROR: Could not find Eigen source directory!"
  exit 1
fi

echo "Configuring Eigen with CMake..."
cmake -S "${EIGEN_SRC}" -B /tmp/eigen-build -DCMAKE_INSTALL_PREFIX="${EIGEN_INSTALL_PREFIX}"

echo "Installing Eigen to ${EIGEN_INSTALL_PREFIX}..."
sudo cmake --install /tmp/eigen-build
echo "Eigen installation complete."

# Optional dependencies:
brew install \
  doxygen \
  ghostscript \
  graphviz
# [installation_documentation]

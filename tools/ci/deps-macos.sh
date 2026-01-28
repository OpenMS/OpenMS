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
  # element.
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
  qt \
  apache-arrow

# Install Eigen 3.4.0 from source (Homebrew's eigen is now 5.x which is incompatible)
# First, remove any Homebrew-installed eigen to avoid version conflicts
echo "Removing any existing Homebrew eigen..."
if command brew list eigen &>/dev/null; then
  echo "Eigen is installed, removing it..."
  command brew uninstall --ignore-dependencies eigen || echo "Warning: Failed to uninstall eigen"
else
  echo "Eigen is not installed via Homebrew, skipping uninstall."
fi

EIGEN_VERSION="3.4.0"
EIGEN_URL="https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz"
echo "Downloading Eigen ${EIGEN_VERSION} from ${EIGEN_URL}..."
curl -fSL --retry 3 --retry-delay 5 "${EIGEN_URL}" -o /tmp/eigen.tar.gz

echo "Extracting Eigen..."
tar -xzf /tmp/eigen.tar.gz -C /tmp

echo "Configuring Eigen with CMake..."
cmake -S /tmp/eigen-${EIGEN_VERSION} -B /tmp/eigen-build -DCMAKE_INSTALL_PREFIX=/opt/homebrew

echo "Installing Eigen..."
sudo cmake --install /tmp/eigen-build
echo "Eigen installation complete."

# Optional dependencies:
brew install \
  doxygen \
  ghostscript \
  graphviz
# [installation_documentation]

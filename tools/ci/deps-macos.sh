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
command brew uninstall --ignore-dependencies eigen 2>/dev/null || true
EIGEN_VERSION="3.4.0"
curl -L "https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz" -o /tmp/eigen.tar.gz
tar -xzf /tmp/eigen.tar.gz -C /tmp
cmake -S /tmp/eigen-${EIGEN_VERSION} -B /tmp/eigen-build -DCMAKE_INSTALL_PREFIX=/opt/homebrew
sudo cmake --install /tmp/eigen-build

# Optional dependencies:
brew install \
  doxygen \
  ghostscript \
  graphviz
# [installation_documentation]

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
  eigen \
  sqlite \
  coinutils \
  cbc \
  cgl \
  clp \
  qtbase \
  qtsvg \
  apache-arrow \
  bash

# Optional dependencies:
brew install \
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

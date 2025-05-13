#!/usr/bin/env bash

set -eu
set -o pipefail

# Code between the following doxygen markers are included in the
# public-facing OpenMS installation instructions.

# [installation_documentation]
# Update the package index:
brew update

# Required dependencies:
brew install --force --overwrite \
  python@3.13 \
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
  qt

# Optional dependencies:
brew install --force --overwrite \
  doxygen \
  ghostscript \
  graphviz
# [installation_documentation]

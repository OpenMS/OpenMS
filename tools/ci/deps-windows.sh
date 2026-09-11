#!/usr/bin/env bash

set -eu
set -o pipefail

# Code between the following doxygen markers are included in the
# public-facing OpenMS installation instructions.

# [installation_documentation]
choco install -y --no-progress \
  cmake \
  ninja

# If you want to install the documentation dependencies:
choco install -y --no-progress \
  graphviz

# Temporary hack to get doxygen installed:
(git clone https://github.com/OpenMS/chocolatey-packages.git &&
  cd chocolatey-packages &&
  make &&
  cd build &&
  choco install doxygen -s .)
# [installation_documentation]

# Use a custom NSIS, which provides 8-k string support and has
# UltraModernUI integrated already:
curl --no-progress-meter -L -o NSIS.tar.gz https://github.com/OpenMS/NSIS/raw/main/NSIS.tar.gz
7z x -so NSIS.tar.gz | 7z x -si -ttar -aoa -o"C:/Program Files (x86)/NSIS/"

# These are only needed in CI:
choco install -y --no-progress \
  ccache \
  rclone

# Use an older version of rsync to match the server:
choco install -y --no-progress --version 6.2.7 rsync

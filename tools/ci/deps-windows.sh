#!/usr/bin/env bash

set -eu
set -o pipefail

choco install -y --no-progress \
  ccache \
  graphviz \
  rsync

# Temporary hack to get doxygen installed:
(git clone https://github.com/OpenMS/chocolatey-packages.git &&
  cd chocolatey-packages &&
  make &&
  cd build &&
  choco install doxygen -s .)

# Use a custom NSIS, which provides 8-k string support and has
# UltraModernUI integrated already:
curl --no-progress-meter -L -o NSIS.tar.gz https://github.com/OpenMS/NSIS/raw/main/NSIS.tar.gz
7z x -so NSIS.tar.gz | 7z x -si -ttar -aoa -o"C:/Program Files (x86)/NSIS/"

#!/usr/bin/env bash

set -eu
set -o pipefail

# Code between the following doxygen markers are included in the
# public-facing OpenMS installation instructions.

# [installation_documentation]
choco install -y --no-progress \
  cmake \
  ninja

# If you want to install the documentation dependencies (Graphviz is only
# needed for the optional 'doc_dot' target, i.e. documentation with all dot
# graphs, and is therefore not installed here: choco install graphviz).
# Temporary hack to get doxygen installed:
(git clone https://github.com/OpenMS/chocolatey-packages.git &&
  cd chocolatey-packages &&
  make &&
  cd build &&
  choco install doxygen -s .)
# [installation_documentation]

# Use a custom NSIS, which provides 8-k string support and has
# UltraModernUI integrated already:
curl --no-progress-meter --fail --retry 3 --retry-all-errors -L -o NSIS.tar.gz https://raw.githubusercontent.com/OpenMS/NSIS/main/NSIS.tar.gz
7z x -so NSIS.tar.gz | 7z x -si -ttar -aoa -o"C:/Program Files (x86)/NSIS/"

# These are only needed in CI. Like doxygen, they come from the packages built
# above (SHA256-checked) rather than from the Chocolatey community feed, which
# fails too often. The rsync package is cwRsync 6.2.7, an older version that
# matches the server:
(cd chocolatey-packages/build &&
  choco install -y --no-progress -s . ccache rclone rsync)

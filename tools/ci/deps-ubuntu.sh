#!/usr/bin/env bash

set -eu
set -o pipefail

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
  --skip-gui-deps)
    # accepted for existing callers; the GUI (and Qt) is built by OpenMS-TOPP-GUI
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
  build-essential \
  cmake \
  git \
  ninja-build \
  patch \
  pkg-config

# Optional dependencies:
sudo apt-get -qq install -y \
  doxygen

# Graphviz is only needed for the optional 'doc_dot' target (documentation
# with all dot graphs), so it is not installed here:
#   sudo apt-get install -y graphviz

# [installation_documentation]

# These are only needed in CI:
sudo apt-get -qq install -y \
  ccache \
  rclone

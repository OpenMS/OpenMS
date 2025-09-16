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
# Update the package index (with error handling for macOS 14):
if ! brew update; then
  # Check if we're on macOS 14 or later where this error is known to occur
  if [[ "$(sw_vers -productVersion | cut -d. -f1)" -ge 14 ]]; then
    echo "Warning: brew update failed on macOS 14+, but continuing with installation..."
    echo "This is a known issue with Homebrew type checking and doesn't affect the build."
  else
    echo "Error: brew update failed. This may indicate a real issue on this macOS version."
    # Exit with error on older macOS versions where update failure is unexpected
    exit 1
  fi
fi

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
  qt \
  apache-arrow

# Optional dependencies:
brew install \
  doxygen \
  ghostscript \
  graphviz
# [installation_documentation]

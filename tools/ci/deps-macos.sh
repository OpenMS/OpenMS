#!/usr/bin/env bash

set -eu
set -o pipefail

# Parse command line arguments
SKIP_DOC_DEPS=false
SKIP_GUI_DEPS=false

while [[ $# -gt 0 ]]; do
  case $1 in
    --skip-doc-deps)
      SKIP_DOC_DEPS=true
      shift
      ;;
    --skip-gui-deps)
      SKIP_GUI_DEPS=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 [--skip-doc-deps] [--skip-gui-deps]"
      exit 1
      ;;
  esac
done

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
  apache-arrow \
  zstd \
  bash

# GUI dependencies (can be skipped for non-GUI builds):
if [ "$SKIP_GUI_DEPS" = false ]; then
  brew install qtsvg
fi

# Optional documentation dependencies:
if [ "$SKIP_DOC_DEPS" = false ]; then
  brew install \
    doxygen \
    ghostscript \
    graphviz
fi
# [installation_documentation]

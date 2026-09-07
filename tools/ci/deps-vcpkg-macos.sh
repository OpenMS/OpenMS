#!/usr/bin/env bash
set -e
set -o pipefail

# Parse command line arguments
SKIP_GUI_DEPS="${SKIP_GUI_DEPS:-false}"
SKIP_DOC_DEPS="${SKIP_DOC_DEPS:-false}"

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

# [installation_documentation]
brew update
brew install \
  ccache \
  ninja \
  bison \
  flex \
  autoconf \
  autoconf-archive \
  automake \
  libtool \
  pkg-config \
  bash \
  uv

# bison/flex: vcpkg requires newer versions than the one Apple ships.
# In CI, prepend the Homebrew versions to PATH via GITHUB_PATH. For 
# local build add equivalent lines:
if [ -n "${GITHUB_PATH:-}" ]; then
  echo "$(brew --prefix bison)/bin" >> "$GITHUB_PATH"
  echo "$(brew --prefix flex)/bin" >> "$GITHUB_PATH"
else
  export PATH="$(brew --prefix bison)/bin:$(brew --prefix flex)/bin:$PATH"
  #To make it permanent, add the same line to ~/.zshrc
fi

# GUI dependencies (can be skipped for non-GUI builds):
if [ "$SKIP_GUI_DEPS" = false ]; then
  brew install qtbase qtsvg
fi

# Optional documentation dependencies:
if [ "$SKIP_DOC_DEPS" = false ]; then
  brew install \
    doxygen \
    graphviz
fi
# [installation_documentation]
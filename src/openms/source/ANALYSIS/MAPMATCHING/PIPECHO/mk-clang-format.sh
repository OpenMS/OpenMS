#!/usr/bin/env bash

# Generates a .clang-format file with some tweaks from the standard
# OpenMS version.
set -eu
set -o pipefail

edits=(
  # More than 80 columns is ridiculous!
  -e 's/^(ColumnLimit:).*$/\1 80/'
)

sed -E "${edits[@]}" \
  "$(git rev-parse --show-toplevel)/.clang-format" \
  >.clang-format

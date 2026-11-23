#!/usr/bin/env bash
# Runs once in the new dev container, from the repository root (see
# devcontainer.json): fetches the submodules the vcpkg presets need and prints
# the preset to use on this machine.

set -euo pipefail

cd "$(dirname "$0")/.."

# A workspace cloned into a container volume, or bind-mounted from a host user
# with another UID, is not owned by the container user; git refuses to work in
# it unless it is marked safe.
git config --global --add safe.directory "$PWD"

# vcpkg: the presets' toolchain file (vcpkg/scripts/buildsystems/vcpkg.cmake).
# THIRDPARTY: the search engines and other tools the TOPP adapters and their
# tests call (shallow, as configured in .gitmodules).
git submodule update --init vcpkg THIRDPARTY

case "$(uname -m)" in
  x86_64|amd64) preset=linux-x64-relwithdebinfo ;;
  aarch64|arm64) preset=linux-arm64-relwithdebinfo ;;
  *) echo "error: no OpenMS preset for machine $(uname -m)" >&2; exit 1 ;;
esac

# Put the THIRDPARTY tools on PATH for new shells. VS Code reads its remote
# environment from an interactive login shell, so CMake Tools sees them too.
marker="# OpenMS THIRDPARTY tools (added by .devcontainer/post-create.sh)"
if ! grep -qxF "$marker" ~/.bashrc 2>/dev/null; then
  cat >>~/.bashrc <<EOF

$marker
for _openms_tp in "$PWD/THIRDPARTY/All"/* "$PWD/THIRDPARTY/Linux/$(uname -m)"/*; do
  [ -d "\$_openms_tp" ] && PATH="\$PATH:\$_openms_tp"
done
unset _openms_tp
export PATH
EOF
fi

cat <<EOF

OpenMS dev container ready. To build (the first configure builds all vcpkg
dependencies, which takes a while; later containers reuse the binary cache in
the openms-vcpkg-cache volume):

  cmake --preset $preset
  cmake --build --preset $preset
  ctest --preset $preset

or select the "$preset" configure preset in CMake Tools.
EOF

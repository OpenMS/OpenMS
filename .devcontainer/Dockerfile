# OpenMS development container: Ubuntu 24.04 with the build tools that the
# vcpkg presets in CMakePresets.json need. vcpkg itself builds the C++
# dependencies (Boost, Arrow, Xerces-C, COIN-OR, ...) on the first configure;
# see .devcontainer/README.md.
#
# The build context is the repository root (see devcontainer.json), restricted
# by Dockerfile.dockerignore to the two setup scripts copied below.
#
# The base image is multi-arch (linux/amd64 and linux/arm64), and both scripts
# pick the architecture from the machine they run on.
FROM mcr.microsoft.com/devcontainers/cpp:1-ubuntu-24.04

COPY tools/ci/deps-ubuntu.sh tools/ci/install_dotnet_sdk_linux.sh /tmp/openms-setup/

# - deps-ubuntu.sh: compilers, CMake, Ninja, autotools (for the coin-or ports),
#   pkg-config, Qt 6 (WITH_GUI is ON by default), Doxygen and ccache.
#   It calls add-apt-repository, which needs software-properties-common.
# - libicu74: needed by the .NET SDK and runtime.
# - install_dotnet_sdk_linux.sh: the pinned, checksum-verified .NET 8 SDK in
#   /usr/share/dotnet. WITH_THERMO_RAW is ON by default, and the Thermo bridge
#   compiles against the SDK's nethost headers.
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
       software-properties-common \
       libicu74 \
    && DEBIAN_FRONTEND=noninteractive bash /tmp/openms-setup/deps-ubuntu.sh \
    && bash /tmp/openms-setup/install_dotnet_sdk_linux.sh /usr/share/dotnet \
    && rm -rf /tmp/openms-setup /var/lib/apt/lists/*

ENV DOTNET_ROOT=/usr/share/dotnet \
    DOTNET_CLI_TELEMETRY_OPTOUT=1 \
    PATH=/usr/share/dotnet:$PATH

# vcpkg's binary cache and download cache. devcontainer.json mounts a named
# volume here, so that rebuilding the container does not rebuild every
# dependency. Docker initialises an empty volume from the image directory,
# including its owner, so the non-root user can write to it.
ENV VCPKG_DEFAULT_BINARY_CACHE=/vcpkg-cache/archives \
    VCPKG_DOWNLOADS=/vcpkg-cache/downloads
RUN mkdir -p /vcpkg-cache/archives /vcpkg-cache/downloads \
    && chown -R vscode:vscode /vcpkg-cache

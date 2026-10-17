# OpenMS dev container

A [dev container](https://containers.dev/) for VS Code, GitHub Codespaces and
other tools that support `devcontainer.json`. It builds OpenMS with the vcpkg
presets from `CMakePresets.json`, as described in
`doc/doxygen/install/install-vcpkg.doxygen`.

## What it contains

- `mcr.microsoft.com/devcontainers/cpp:1-ubuntu-24.04` (linux/amd64 and
  linux/arm64), plus the packages from `tools/ci/deps-ubuntu.sh`: compilers,
  CMake, Ninja, autotools, Qt 6, Doxygen and ccache.
- The .NET 8 SDK in `/usr/share/dotnet`, from
  `tools/ci/install_dotnet_sdk_linux.sh`. `WITH_THERMO_RAW` is ON by default,
  and the Thermo RAW bridge compiles against the SDK's nethost headers.
- A named volume, `openms-vcpkg-cache`, at `/vcpkg-cache`, which holds vcpkg's
  binary cache and downloads. The first configure builds every dependency
  (this takes a while); rebuilt containers reuse the cache. Remove the volume
  with `docker volume rm openms-vcpkg-cache` to start over.

On creation, `post-create.sh` initializes the `vcpkg` and `THIRDPARTY`
submodules and puts the THIRDPARTY tools on `PATH`.

## Building

The default presets build RelWithDebInfo, which can be debugged and runs the
tests at a reasonable speed:

| Machine                            | Preset                       |
|------------------------------------|------------------------------|
| x86_64 (most PCs, Codespaces)      | `linux-x64-relwithdebinfo`   |
| arm64 (Apple silicon, ARM servers) | `linux-arm64-relwithdebinfo` |

```bash
cmake --preset linux-x64-relwithdebinfo
cmake --build --preset linux-x64-relwithdebinfo
ctest --preset linux-x64-relwithdebinfo
```

CMake Tools is set to use presets: pick the configure preset in its status bar
or with *CMake: Select Configure Preset*. The `-debug` and `-release` presets
work as well and share the vcpkg binary cache, since the Linux triplets build
both configurations of each dependency. The build directory is
`build/<preset>`.

The container needs at least 4 CPUs, 16 GB of memory and 64 GB of disk; give
Docker more CPUs to build faster.

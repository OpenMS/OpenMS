GNU/Linux
=========================


```{include} installation-with-conda.md
:start-after: "% start-after"
```

## Install via package managers

Packaged versions of **OpenMS** are provided for Fedora, OpenSUSE, Debian, and Ubuntu. You can find them to download
[here](https://pkgs.org/download/openms). For other GNU/Linux distributions or to obtain the most recent version of the
library, installation should be done via building from the source code.

```{important}
These packages are not directly maintained by the OpenMS team and they can not be guaranteed to have the
same behaviour as when building it from source code. Also, their availability and version is subject to change and
support might be limited (due to unforeseen or untested behaviour). It is suggested not to install them parallel to our
Debian package.
```

```{note}
Some thirdparty software used via adapter tools in OpenMS might also require an installed JavaVM.
```

## Install via the provided Debian package

For Debian-based Linux users, it is suggested to use the Debian package attached to each
[OpenMS release](https://github.com/OpenMS/OpenMS/releases/latest), built for x86_64 and for aarch64 (ARM64).
It is tested on Ubuntu 24.04 and needs glibc 2.38 or newer, so it does not install on Ubuntu 22.04 or Debian 12. On
older distributions, use conda (see above) or a container (see below).

Install it with `apt`, which resolves its dependencies, such as Qt 6, from your distribution's repositories:

```bash
sudo apt install ./OpenMS-<version>-Debian-Linux-<architecture>.deb
```
If you encounter errors, troubleshoot using the following steps.

1. Packages are missing, for example the Qt 6 libraries.

   The package depends on the libraries of current distributions. Older distributions do not provide them; see the
   requirements above.

2. Error while executing a tool

   To ensure the tool functionality, make sure you add the `OPENMS_DATA_PATH` variable to your environment as follow
   `export OPENMS_DATA_PATH=/usr/share/OpenMS`

3. Adapters are not finding thirdparty applications

   Executables for thirdparty applications can be found in:
   `/usr/share/OpenMS/THIRDPARTY`
   Add the folders in your `PATH` for a convenient use of the adapters.

```{include} run-in-container.md
:start-after: "% start-after"
```

## Reading Thermo Fisher RAW files

OpenMS reads Thermo Fisher `.raw` files natively through the openms-thermo-bridge, which is
enabled by default in the release binaries on supported platforms. This requires a **.NET 8
runtime** to be present at run time so that the managed bridge libraries can be loaded.

Install it from the [.NET download page](https://dotnet.microsoft.com/download) or via your
distribution's package manager, for example:

```bash
# Debian/Ubuntu
sudo apt-get install dotnet-runtime-8.0
```

If .NET is installed to a non-standard location (for example via the `dotnet-install.sh`
script), point the `DOTNET_ROOT` environment variable at the install directory — the folder
that contains the `dotnet` host and the `shared/` sub-directory — so the bridge can locate the
runtime:

```bash
export DOTNET_ROOT=/usr/share/dotnet
```

## Build OpenMS from source

To build OpenMS from source, follow the build instructions for [Linux](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/install_linux.html).

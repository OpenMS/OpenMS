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

For Debian-based Linux users, it is suggested to  use the [deb-package](https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/latest/) provided. It is most easily installed with **[gdebi](https://launchpad.net/gdebi)**
which automatically resolves the dependencies available in the PPA Repositories.

```bash
sudo apt-get install gdebi
sudo gdebi /PATH/TO/OpenMS.deb
```
If you encounter errors with unavailable packages, troubleshoot using the following steps.

1. Qt5 (or one of its packages, e.g. `qt5xbase`) is missing.

   It might be because your Debian is too old to have a recent enough version in its official repositories. It is
   suggested to use the same packages that are used while building (make sure to adapt the Qt version and your
   Debian/Ubuntu version, here Xenial):
   ```bash
   sudo add-apt-repository ppa:beineri/opt-qt59-xenial
   sudo apt-get update
   ```
   Run the installation again.

2. ICU with its `libicu` is missing.

   You can find the missing version on [pkgs.org](https://pkgs.org) and install it with `gdebi`, too. You can have
   multiple versions of ICU installed.

3. Error while executing a tool

   To ensure the tool functionality, make sure you add the `OPENMS_DATA_PATH` variable to your environment as follow
   `export OPENMS_DATA_PATH=/usr/share/OpenMS`

4. Thirdparty installation of Qt5 in step 1

   Make sure you source the provided environment file using:
   `source /opt/qt59/bin/qt59-env.sh`

5. Adapters are not finding thirdparty applications

   Executables for thirdparty applications can be found in:
   `/usr/share/OpenMS/THIRDPARTY`
   Add the folders in your `PATH` for a convenient use of the adapters.

```{include} run-in-container.md
:start-after: "% start-after"
```

## Reading Thermo Fisher RAW files

OpenMS reads Thermo Fisher `.raw` files natively through the openms-thermo-bridge, which is
enabled by default in the release binaries on supported platforms. This requires a **.NET 8
(or newer) runtime** to be present at run time so that the managed bridge libraries can be
loaded.

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

```{note}
Native Thermo RAW reading is not available on Linux/aarch64 (ARM64), because Thermo's
RawFileReader does not ship native dependencies for that platform.
```

## Build OpenMS from source

To build OpenMS from source, follow the build instructions for [Linux](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/install_linux.html).

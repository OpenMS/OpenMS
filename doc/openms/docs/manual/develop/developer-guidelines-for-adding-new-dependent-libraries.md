Developer Guidelines For Adding New Dependent Libraries
======================================================

## Our dependency library philosophy

In short, requirements for adding a new library are:
- indispensable functionality
- license compatibility
- availability for all platforms

## Indispensable functionality

In general, adding a new dependency library (of which we currently have more than a handful, e.g. Xerces-C or ZLib)
imposes a significant integration and maintenance effort. Thus, the new library should add 
**indispensable functionality**. If the added value does not compensate for the overhead, alternative solutions 
encompass:


- write it yourself and add to the OpenMS library (i.e. its repository) directly
- write a TOPPAdapter which calls an external executable (placing the burden on the user to supply the executable)

## License compatibility

OpenMS has a BSD-3 clause license and we try hard to remove dependencies of GSL-like libraries. Therefore, a new library
with e.g. LGPL-2 would be prohibitive.


## C++ standard compatibility

New dependency libraries needs to be compatible and therefore compilable with the same C++ standard as OpenMS.

## Availability for all platforms

OpenMS has been designed for Windows, macOS, and Linux. Therefore, the new dependency library needs to be designed for
these platforms.

OpenMS obtains its dependencies through [vcpkg](https://vcpkg.io) on all three platforms: the presets in
`CMakePresets.json` set `OPENMS_USE_VCPKG=ON`, and vcpkg builds the libraries listed in the manifest `vcpkg.json`.
The library therefore needs a vcpkg port that builds on every triplet OpenMS uses (see `CMakePresets.json` and
`vcpkg-overlays/triplets`), or you have to write one (see below).

- on **Windows** OpenMS uses the `x64-windows-static-md` triplets: the library is linked statically, against the
  **dynamic** VS-C++ runtime, like OpenMS itself. The debug preset builds every dependency in a debug and a release
  variant, so the debug and release runtimes are never mixed; the port does not need to take care of that itself.
  The library must build with the minimum Visual Studio version OpenMS supports.

- on **macOS** it should be ensured that the library can be built on recent macOS versions with Apple Clang and the
  mac specific _libc++_. Ideally the package is also available via **HomeBrew**, which the macOS pyOpenMS wheel and
  builds without vcpkg use.

- on **Linux** since we (among other distributions) feature an OpenMS Debian package, and the images in `dockerfiles/`
  build OpenMS without vcpkg from distribution packages, the new library should be available as Debian/Ubuntu package
  as well, or be linked statically during the OpenMS packaging build.

## How to add it to the build

1. **Add it to the vcpkg manifest.** Add the port to `vcpkg.json`. If the dependency is optional, add it as a feature
   under `"features"` instead of to the top-level `"dependencies"`, and pair it with a CMake option `WITH_<NAME>`
   (see e.g. the `hdf5` feature and `WITH_HDF5`). Enabling such a dependency needs both
   `-DWITH_<NAME>=ON` and `-DVCPKG_MANIFEST_FEATURES="...;<feature>"`. If CI should build with it, add the feature to
   `VCPKG_MANIFEST_FEATURES` of the `*-ci` presets in `CMakePresets.json`.

2. **Add an overlay port if needed.** If the library is not in the vcpkg registry, or the registry's port needs patches
   or different build options, add a port under `vcpkg-overlays/ports/<name>/` (a `vcpkg.json`, a `portfile.cmake`, and
   patches created with `git diff` or `diff -Naur`). A port there replaces the registry port of the same name.
   `vcpkg-configuration.json` already points vcpkg to this directory. See the
   [vcpkg build guide](https://archive.openms.de/openms/Documentation/nightly/latest/html/install_vcpkg.html)
   and the vcpkg documentation on [overlay ports](https://learn.microsoft.com/en-us/vcpkg/concepts/overlay-ports).

3. **Find it in CMake.** Add the `find_package()` call and the target to link to `cmake/cmake_findExternalLibs.cmake`,
   guarded by the `WITH_<NAME>` option for an optional dependency. Prefer the library's own CMake config package
   (`find_package(<Name> CONFIG)`); add a find module under `cmake/Modules/` only if the library has none.

4. **Keep builds without vcpkg working.** The same `find_package()` call must also find the library when OpenMS is
   configured with `-DOPENMS_USE_VCPKG=OFF` against system packages (apt, Homebrew, conda). Add the package to
   `tools/ci/deps-*.sh` where applicable and to the images in `dockerfiles/`.

Then test the build on your platform, including a fresh configure with a preset so that vcpkg builds the new port.
Make sure the library is correctly shipped in the installer packages and in the pyOpenMS wheels (especially shared
libraries, and especially on Windows), and that its license is included.

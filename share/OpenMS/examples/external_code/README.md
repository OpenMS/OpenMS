# Example CMake project using OpenMS

Example project for external code using OpenMS library and headers.
You can modify the build system via CMakeLists.txt, e.g., to add more C++ classes, alter build flags, or add additional dependencies.

`find_package(OpenMS CONFIG)` provides the OpenMS library as the imported target `OpenMS::OpenMS` (the OpenSWATH
algorithm library as `OpenMS::OpenSwathAlgo`, and the TOPP tool framework as `OpenMS::OpenMS_CLI`, which a program
deriving from `TOPPBase` links instead; it is part of an installation that includes the CLI layer, so such a
program requests it with `COMPONENTS CLI`); linking it supplies include directories, compile features and
dependencies. Building OpenMS and consuming its CMake package both require CMake 3.24 or newer.

## Usage

Assuming everything happens in the directory `~/dev`, e.g., the OpenMS sources are located in `~/dev/OpenMS` and OpenMS was compiled in `~/dev/OpenMS-build`. For any details on how to compile OpenMS please check out either the documentation shipped with OpenMS or online at http://www.openms.de/documentation.

 1. Compile OpenMS (e.g., in `~/dev/OpenMS-build`)
 2. Create a new build directory for this example project (e.g., `mkdir ~/dev/example-build`)
 3. Call CMake from within the new directory, with the source dir for the external code as last argument (e.g., `cmake -G "<generator used for OpenMS>" ~/dev/OpenMS/share/OpenMS/examples/external_code/`). You can also copy `~/dev/OpenMS/share/OpenMS/examples/external_code/` to any other place and reference that instead.
    If OpenMS was built with vcpkg (e.g. with `cmake --preset <preset>`), pass the same vcpkg toolchain and installation, so that CMake finds the libraries OpenMS depends on (e.g., Boost):
    `-DCMAKE_TOOLCHAIN_FILE=$HOME/dev/OpenMS/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_INSTALLED_DIR=$HOME/dev/OpenMS/build/<preset>/vcpkg_installed -DVCPKG_TARGET_TRIPLET=<triplet of the preset>`.
 
**Note**: In general you should try to use the same setup (compiler etc.) for OpenMS and your project. Especially on Windows you need to use the same CMake Generator for OpenMS and the new project. 
 
Should the above step fail because CMake couldn't find OpenMS you can specify the OpenMS instance when calling CMake, e.g., 

```cmake -G "<generator used for OpenMS>" -D OpenMS_DIR=~/Development/OpenMS-build/ ~/Development/OpenMS/share/OpenMS/examples/external_code/```

For an installed OpenMS, `OpenMS_DIR` is the directory holding `OpenMSConfig.cmake`, i.e. `<prefix>/lib/cmake/OpenMS`
(`<prefix>/CMake` on Windows); alternatively add `<prefix>` to `CMAKE_PREFIX_PATH`.

# External project example

This project shows how to compile custom code against an installed OpenMS: `find_package(OpenMS CONFIG)`
provides the imported targets `OpenMS::OpenMS`, `OpenMS::OpenSwathAlgo` and `OpenMS::OpenMS_CLI` (the TOPP tool
framework; `TestExternalCodeCLI` derives a tool from `TOPPBase` against it), and `OpenMS::OpenMS_GUI` for
an installation built with `WITH_GUI=ON`. Requesting the `GUI` component (`COMPONENTS GUI` or
`OPTIONAL_COMPONENTS GUI`) additionally finds the Qt6 modules the GUI library links and sets
`OpenMS_GUI_FOUND`; a project that links `OpenMS::OpenMS_GUI` without requesting the component has to
find those Qt6 modules itself. The un-namespaced names `OpenMS`, `OpenSwathAlgo`, `OpenMS_CLI` and `OpenMS_GUI` of
earlier releases remain available as aliases. Consuming projects need CMake 3.22 or newer.

It also serves as the test that the CMake package of an OpenMS installation works: when OpenMS is
configured with `-DOPENMS_TEST_INSTALLED_CONSUMER=ON` (on in the CI presets), the CTest tests
`TestExternalCode_*` install the development components into `<build>/installed-consumer/prefix`,
then configure, build and run this project against that installation.

## Usage

 1. Build and install OpenMS (e.g. `cmake --install OpenMS-build --prefix ~/OpenMS-install`).
 2. Create a build directory for this project (e.g. `mkdir ~/example-build`).
 3. Configure it against the installation, using the same compiler, generator and dependency
    provider (vcpkg toolchain file or `CMAKE_PREFIX_PATH`) as OpenMS:

        cmake -S <OpenMS>/src/tests/external -B ~/example-build -G "<generator used for OpenMS>" \
              -DOpenMS_DIR="$HOME/OpenMS-install/lib/cmake/OpenMS"

    Alternatively add the installation prefix to `CMAKE_PREFIX_PATH` instead of setting
    `OpenMS_DIR`. On Windows the CMake package lives in `<prefix>/CMake`.
 4. Build and run the tests: `cmake --build ~/example-build` and `ctest --test-dir ~/example-build`.

On Windows, `<prefix>/bin` (the OpenMS DLLs) and the directories of the dependency DLLs must be
on `PATH` when running your executables; this project's `CMakeLists.txt` shows how to do that
for CTest using `OPENMS_LIB_DIR` from `OpenMSConfig.cmake`.

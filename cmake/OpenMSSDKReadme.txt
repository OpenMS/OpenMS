OpenMS SDK
==========

This archive holds what you need to build your own C++ code against OpenMS:

  include/            headers of OpenMS, OpenSwathAlgo, the TOPP tool framework
                      (OpenMS_CLI) and the bundled third-party headers
  lib/ (bin/ on Windows)
                      the OpenMS libraries and the shared libraries they depend on
  lib/cmake/OpenMS/ (CMake/ on Windows)
                      the CMake package (OpenMSConfig.cmake)
  share/OpenMS/       shared data (CV files, chemistry, databases) needed at runtime;
                      share/OpenMS/models holds the PeptDeep ONNX models for the
                      PeptDeep inference classes, PeptDeepRTInference, PeptDeepCCSInference
                      and PeptDeepMS2Inference (the ONNX Runtime
                      library is bundled with the other dependencies)

It is relocatable: extract it anywhere.


Using it from CMake (3.24 or newer)
-----------------------------------

  find_package(OpenMS CONFIG REQUIRED)              # OpenMS::OpenMS, OpenMS::OpenSwathAlgo
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
                                                    # + OpenMS::OpenMS_CLI (TOPPBase-style tools)

  A version request has to name major and minor version (find_package(OpenMS @OPENMS_VERSION_MAJOR_MINOR@ ...)):
  the package is compatible within one minor version only.

  add_executable(mytool mytool.cpp)
  target_link_libraries(mytool PRIVATE OpenMS::OpenMS_CLI)

Configure your project with the SDK on the prefix path:

  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/OpenMS-SDK-<version>-<platform>

A complete example project (library, plain program and TOPP-style tool) is in
the OpenMS sources under src/tests/external.


Requirements of your project
----------------------------

  * @BOOST_REQUIREMENT@
  * No Qt: the SDK holds the core and CLI layers, which do not use it. The GUI
    library (OpenMS::OpenMS_GUI) is not part of the SDK; build OpenMS from source
    if you need it.
  * The same compiler family the SDK was built with:
      - Windows: MSVC (Visual Studio 2022 17.14 or newer), x64, Release configuration
        with the dynamic runtime (/MD). A Debug build of your code (/MDd) must not
        be linked against this SDK; build OpenMS from source for that.
      - macOS: Apple Clang, arm64.
      - Linux: GCC with a glibc and libstdc++ at least as new as those of
        Ubuntu 24.04, on which the SDK is built.


Running your programs
---------------------

  * Linux and macOS: the libraries find each other and their bundled
    dependencies by themselves; CMake gives your executables an RPATH to the
    SDK's lib/ directory in the build tree.
  * Windows: put the SDK's bin/ directory on PATH.
  * Set OPENMS_DATA_PATH to the SDK's share/OpenMS directory: OpenMS looks for
    its shared data next to the running executable, which for your programs is
    not inside the SDK.

Other ways to get OpenMS for development: the conda package libopenms
(bioconda, Linux and macOS) or building OpenMS from source
(https://openms.readthedocs.io).

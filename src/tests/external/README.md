# External project example

This project shows how to compile custom code against OpenMS.

Additionally this project serves as a test case for nightly builds to ensure that the prorcess of finding the CMake exports of OpenMS works and compiling against them works.

## Usage

Assuming everything happens in the directory `~/Development`, e.g., the OpenMS sources are located in `~/OpenMS` and OpenMS was compiled in `~/OpenMS-build`. For any details on how to compile OpenMS please check out either the documentation shipped with OpenMS or online at http://www.openms.de/documentation.

 1. Compile OpenMS (e.g., in `~/Development/OpenMS-build`)
 2. Create the environment variable `OPENMS_BUILD_TREE` that points to the build directory OpenMS (e.g., `export OPENMS_BUILD_TREE=~/Development/OpenMS-build` or `set OPENMS_BUILD_TREE=~/Development/OpenMS-build`). **Note**: This is not necessary in general but is used in the nightly build/test setting to make sure that we use a specific OpenMS instance.
 3. Create a new build directory for this example project (e.g., `mkdir ~/Development/example-build`)
 4. Call CMake from the new directory and add the path to this directory (e.g., `cmake -G "<generator used for OpenMS>" ~/Development/OpenMS/share/OpenMS/examples/external_code/`). 
 
**Note**: In general you should try to use the same setup (compiler etc.) for OpenMS and your project. Especially on Windows you need to use the same CMake Generator for OpenMS and the new project.

# Example CMake project using OpenMS

Example project for external code using OpenMS library and headers.
You can modify the build system via CMakeLists.txt, e.g., to add more C++ classes, alter build flags, or add additional dependencies.

## Usage

Assuming everything happens in the directory `~/dev`, e.g., the OpenMS sources are located in `~/dev/OpenMS` and OpenMS was compiled in `~/dev/OpenMS-build`. For any details on how to compile OpenMS please check out either the documentation shipped with OpenMS or online at http://www.openms.de/documentation.

 1. Compile OpenMS (e.g., in `~/dev/OpenMS-build`)
 2. Create a new build directory for this example project (e.g., `mkdir ~/dev/example-build`)
 3. Call CMake from within the new directory, referencing where the OpenMS contrib was build (for external dependencies like xerces), and as last argument the source dir for the external code (e.g., `cmake -G "<generator used for OpenMS>" -DOPENMS_CONTRIB_LIBS=/path/to/openms_contrib_build  ~/dev/OpenMS/share/OpenMS/examples/external_code/`). You can also copy `~/dev/OpenMS/share/OpenMS/examples/external_code/` to any other place and reference that instead.
 
**Note**: In general you should try to use the same setup (compiler etc.) for OpenMS and your project. Especially on Windows you need to use the same CMake Generator for OpenMS and the new project. 
 
Should the above step fail because CMake couldn't find OpenMS you can specify the OpenMS instance when calling CMake, e.g., 

```cmake -G "<generator used for OpenMS>" -D OpenMS_DIR=~/Development/OpenMS-build/ ~/Development/OpenMS/share/OpenMS/examples/external_code/```

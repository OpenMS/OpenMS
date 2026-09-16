External Code using OpenMS
==========================

If OpenMS' TOPP tools are not enough in a certain scenario, you can either request a change to OpenMS, if you
feel this functionality is useful for others as well, or modify/extend OpenMS privately. For the latter, there are 
multiple ways to do this:

- Modify the developer version of OpenMS by changing existing tools or adding new ones.
- Use an **External Project** to write a new tool, while not touching OpenMS itself (see below on how to do that).

Once you've finished your new tool, and it only needs to run on the development machine. To ship it to a new client 
machine, see, read further in this document.

## Compiling external code

It is very easy to set up an environment to write your own programs using OpenMS. Build OpenMS from source or
install a development package, then let CMake locate it with `find_package(OpenMS CONFIG)`. The package provides
the OpenMS libraries as imported targets, which carry the include directories, compile features and dependencies
your code needs:

- `OpenMS::OpenMS`: the OpenMS library (link this)
- `OpenMS::OpenSwathAlgo`: the OpenSWATH algorithm library
- `OpenMS::OpenMS_CLI`: the TOPP tool framework (`TOPPBase`, `ToolHandler`, `INIUpdater`, ...); link this instead of
  `OpenMS::OpenMS` when your program derives from `TOPPBase` (it carries `OpenMS::OpenMS` transitively), and request
  the `CLI` component (`find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)`) so that an installation without the
  framework is rejected when your project is configured
- `OpenMS::OpenMS_GUI`: the GUI library of an OpenMS built with `WITH_GUI=ON` whose installation includes it;
  request the `GUI` component (`find_package(OpenMS CONFIG REQUIRED COMPONENTS GUI)`) to also find the Qt6 modules
  it was built against, otherwise a project linking it has to find those Qt6 modules itself

The installed package is layered. The core layer (the install components `library`, `OpenMS_headers`,
`OpenSwathAlgo_headers`, `thirdparty_headers`, `share` and `cmake`) is always present; the TOPP tool framework
(`library_cli`, `cmake_cli` and `OpenMS_CLI_headers`) and the GUI library (`library_gui`, `cmake_gui` and
`OpenMS_GUI_headers`) are optional layers on top of it, each linking the layer below. `cmake --install <build>`
installs everything, `cmake --install <build> --component <name>` one component at a time; the pyOpenMS wheels, for
instance, are built against a core-only installation. `OpenMSConfig.cmake` provides the targets of the layers an
installation contains and reports them with `OpenMS_CLI_FOUND` and `OpenMS_WITH_GUI`.

The un-namespaced names `OpenMS` and `OpenSwathAlgo` remain available as aliases for projects written against earlier
releases, as do `OpenMS_CLI` and `OpenMS_GUI` when their layers are installed. The package also reports the version of the installation (`OpenMS_VERSION`), its layers
(`OpenMS_CLI_FOUND`, `OpenMS_WITH_GUI`), its build options (`OpenMS_WITH_HDF5`, `OpenMS_WITH_OPENTIMS`,
`OpenMS_WITH_THERMO_RAW`, `OpenMS_WITH_OPENMP`, `OpenMS_BUILD_TOPP_TOOLS`) and its directories (`OPENMS_DATA_DIR`,
`OPENMS_LIB_DIR`, `OPENMS_BIN_DIR`, `OPENMS_DOC_DIR`). Consuming projects need CMake 3.19 or newer.

```{note}
CMake finds OpenMS through `OpenMS_DIR`, the directory holding `OpenMSConfig.cmake`: `<prefix>/lib/cmake/OpenMS` of an
installation (`<prefix>/CMake` on Windows) or the OpenMS build directory. Alternatively add `<prefix>` to
`CMAKE_PREFIX_PATH`. Use the same compiler, generator and dependency locations (e.g. `OPENMS_CONTRIB_LIBS`) as for
the OpenMS build.
```

The example that follows will be explained in details:

```cmake
cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

### example CMakeLists.txt to develop C++ programs using OpenMS
project("Example_Project_using_OpenMS")

## list all your executables here (a corresponding .cpp file should exist, e.g. Main.cpp)
set(my_executables
  Main
)

## list all classes here, which are required by your executables
## (all these classes will be linked into a library)
set(my_sources
  ExampleLibraryFile.cpp
)

## find OpenMS: provides the imported targets OpenMS::OpenMS (the library) and OpenMS::OpenSwathAlgo, plus
## OpenMS::OpenMS_CLI (the TOPP tool framework, for programs deriving from TOPPBase) when the installation
## includes it (add COMPONENTS CLI to require it).
## If this fails, point CMake at an installation or build tree with -DOpenMS_DIR=<prefix>/lib/cmake/OpenMS
## (<prefix>/CMake on Windows, or the OpenMS build directory), or add <prefix> to CMAKE_PREFIX_PATH.
find_package(OpenMS CONFIG REQUIRED)
message(STATUS "Found OpenMS ${OpenMS_VERSION} at ${OpenMS_DIR}")

## library with additional classes from above
add_library(my_custom_lib STATIC ${my_sources})
target_link_libraries(my_custom_lib PUBLIC OpenMS::OpenMS)

## add targets for the executables
foreach(i ${my_executables})
  add_executable(${i} ${i}.cpp)
  ## link executables against OpenMS
  target_link_libraries(${i} PRIVATE OpenMS::OpenMS my_custom_lib)
endforeach(i)
```

The command `project` defines the name of the project, the name is only of interest of you're working in an IDE or want
to export this project's targets. To compile the program, append it to the `my_executables` list. If you use object
files (classes which do not contain a main program), append them to the `my_sources` list. In the next step CMake
creates a statically linked library of the object files, listed in `my_sources`, and links it against `OpenMS::OpenMS`;
the executables inherit that dependency. This simple CMakeLists.txt example can be extended to also build shared
libraries, include other external libraries and so on.

An example external project can be found in `OpenMS/share/OpenMS/examples/external_code`. Copy these files to a separate
directory and use CMake to configure it (here as an in-source build).

```bash
cd <path_to_external_project>
cmake -G "<generator>" .
```

For more information visit the website of cmake at cmake.org and consult the documentation.

```{important}
Have fun coding with OpenMS!
```

## Shipping external code to a new machine

If you've modified OpenMS itself and not used an external project use our installer scripts, to build your own OpenMS
installer for your platform (see our internal FAQ which is built using "make doc_internal") and ship that to a client
machine.

If you've used an external project and have a new executable (+ an optional new library), use the installer approach as
well, and manually copy the new executable to the `TOPP` binary directory (e.g. on Windows this could be
`c:/program files/OpenMS/bin`, on Linux it could be `/bin`.

If you do not use the installer, copy all required files manually, plus a few extra steps, see below. What needs to be
done is a little platform dependent, thus very cumbersome to explain. Look at the cmake installer scripts, to see whats
required (for macOS and Linux see `OpenMS/cmake/package*.cmake`).

In short:

- copy the `OpenMS/share/OpenMS` directory to the client machine (e.g `<client/my_dir>/share`) and set the environment
  variable `OPENMS_DATA_PATH` to this directory
- copy the OpenMS library (`OpenMS.dll` for Windows or `OpenMS.so/.dylib` for Linux/macOS) to `<client/my_dir>/bin`.
- copy all Qt4 libraries to the client `<client/my_dir>/bin` or on Linux/macOS make sure you have installed the Qt4 
  package.
- [Windows only] copy Xerces dll (see `contrib/lib`) to `<client/my_dir>/bin`
- [Windows only] install the VS redistributable package (see Microsoft Homepage) on the client machine which corresponds
  to the VS version that was used to compile your code (use the correct redistributable package!, i.e., architecture
  32|64bit, VS version, VS Service Pack version). If you choose the wrong redistributable package, you will get
  "Application failed to initialize properly..." error messages.

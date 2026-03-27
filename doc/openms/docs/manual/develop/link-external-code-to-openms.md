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

It is very easy to set up an environment to write your own programs using OpenMS. Make sure to downloaded and installed
the source package of OpenMS/TOPP properly.

```{note}
You cannot use the `install` target when working with the development version of OpenMS, it must be built and used 
within the build tree.
```

All important compiler settings and preprocessor definitions along with the OpenMS library are available. The most
important variables are:

- `OpenMS_INCLUDE_DIRECTORIES`:  all include directories containing OpenMS headers
- `OPENMS_ADDCXX_FLAGS`: preprocessor macros we require written as `(-DMACRO1 -DMACRO2)`

and the OpenMS target itself (which you can link against).

The example that follows will be explained in details:

```cpp
### example CMakeLists.txt to develop C++ programs using OpenMS
project("Example_Project_using_OpenMS")
cmake_minimum_required(VERSION 3.0)

## list all your executables here (a corresponding .cpp file should exist, e.g. Main.cpp)
set(my_executables
  Main
)

## list all classes here, which are required by your executables
## (all these classes will be linked into a library)
set(my_sources
  ExampleLibraryFile.cpp
)

## find OpenMS configuration and register target "OpenMS" (our library)
find_package(OpenMS)
## if the above fails you can try calling cmake with -D OpenMS_DIR=/path/to/OpenMS/
## or modify the find_package() call accordingly
## find_package(OpenMS PATHS "</path/to/OpenMS//")

# check whether the OpenMS package was found
if (OpenMS_FOUND)
  message(STATUS "\nFound OpenMS at ${OpenMS_DIR}\n")

  ## library with additional classes from above
  add_library(my_custom_lib STATIC ${my_sources})

  ## add targets for the executables
  foreach(i ${my_executables})
    add_executable(${i} ${i}.cpp)
    ## link executables against OpenMS
    target_link_libraries(${i} OpenMS my_custom_lib)
  endforeach(i)


else(OpenMS_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake file not found!")
endif(OpenMS_FOUND)
```

The command `project` defines the name of the project, the name is only of interest of you're working in an IDE or want
to export this project's targets. To compile the program, append it to the `my_executables` list. If you use object 
files (classes which do not contain a main program), append them to the `my_sources` list. In the next step CMake 
creates a statically linked library of the object files, listed in `my_sources`. This simple CMakeLists.txt example can 
be extended to also build shared libraries, include other external libraries and so on.

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

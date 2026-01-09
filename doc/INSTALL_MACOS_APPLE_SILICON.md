# Building OpenMS on macOS(with apple silicon)

This document describes the process that has been used to set up and build OpenMS

On macOS with apple Silicon(arm64) using homebrew. It is targeted towards general issues that new commiters might face: may face.

# Tested environment

- macOS: Apple Silicon (arm64)
- Compiler: AppleClang (Xcode Command Line Tools)
- Package manager: Homebrew
- RAM: 8 GB
- Disk space available during build: ~25 GB

# Dependencies (homebrew)

The following dependencies were installed via Homebrew:

cmake  
qt@6  
boost  
eigen  
xerces-c  
glpk  
zstd  
icu4c

# CMake Configuration

A minimal config that was working on macOS on Apple Silicon:

mkdir build
cd build

cmake .. \
  -DMT_ENABLE_OPENMP=OFF \
  -DQt6_DIR=$(brew --prefix qt@6)/lib/cmake/Qt6 \
  -DCMAKE_PREFIX_PATH="$(brew --prefix qt@6);$(brew --prefix zstd);$(brew --prefix icu4c)"

OpenMP is turned off since the Clang compiler that comes with macOS doesn’t support it out of the box — you’d need to install extra libraries to enable it.

# Disk Space Considerations

A full build of OpenMS on a Mac may consume quite a lot of disk space because of object files and temporary compilation results.

With a free disk space of approximately 25 GB, a build may fail with errors such as:

- Unable to rename temporary object files
- Missing intermediate CMake files

A minimum of 40 GB of free space is required for a stable installation.

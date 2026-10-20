# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# WebAssembly build of the OpenMS libraries with Emscripten.
#
# Included by the top-level CMakeLists.txt right after project() when the
# Emscripten toolchain is in use (CMake preset wasm32-emscripten-release), before
# the build options are declared, so that the defaults set here are the ones the
# option() calls find in the cache. See src/openms_wasm/README.md.
#
# A WebAssembly module is a 32-bit (wasm32) program without processes, sockets or
# shared libraries, and, unless built with OPENMS_WASM_THREADS, without threads.
# Running in a browser it also has no access to the file system of the machine:
# files live in Emscripten's virtual file system. Everything that needs one of
# those is switched off; what remains is the OpenMS library (OpenSwathAlgo,
# OpenMS, OpenMS_CLI), built as static libraries, and the JavaScript bindings of
# src/openms_wasm.
#------------------------------------------------------------------------------

message(STATUS "Emscripten ${EMSCRIPTEN_VERSION}: building OpenMS for WebAssembly (${CMAKE_SYSTEM_PROCESSOR})")

## Off by default; they can be enabled, but are not part of a library build.
##  - TOPP tools and class tests become Node.js programs that use the host file
##    system (see openms_add_executable_compiler_flags in compiler_flags.cmake).
##  - The documentation needs doxygen and runs the (native) TOPP tools.
set(BUILD_TOPP_TOOLS OFF CACHE BOOL "Build the TOPP command-line applications.")
set(ENABLE_CLASS_TESTING OFF CACHE BOOL "Enables tests for library classes. Should be disabled only on time constraints (e.g. chunking during continuous integration).")
set(ENABLE_DOCS OFF CACHE BOOL "Indicates whether documentation should be built.")
set(HAS_XSERVER OFF CACHE BOOL "Indicates if an X server is available. If set to Off it will disable certain tests and the doc target.")
set(INSTALL_OPENMS_EXAMPLES OFF CACHE BOOL "Install the OpenMS example data.")

## The LP solver: COIN-OR needs a (Fortran) LAPACK, which is not built for WebAssembly;
## the vcpkg manifest provides HiGHS there instead.
set(LP_SOLVER "HIGHS" CACHE STRING "LP solver to use: AUTO, COIN, GLPK, or HIGHS")

## The bundled SQLite (src/openms/extern/SQLiteCpp) loads extensions with dlopen(),
## which a WebAssembly module without dynamic linking does not have.
set(SQLITE_OMIT_LOAD_EXTENSION ON CACHE BOOL "Enable omit load extension")

## Threads (CMake preset wasm32-emscripten-mt-release): OpenMP, with Emscripten's port
## of LLVM's libomp, and std::thread run on WebAssembly threads, i.e. Web Workers (Node.js:
## worker threads) sharing the memory of the module. Every object of such a module has
## to be compiled with -pthread, the dependencies included: they come from the vcpkg
## triplet wasm32-emscripten-mt. In browsers, shared memory needs a cross-origin isolated
## page; the JavaScript module runs single-threaded elsewhere (src/openms_wasm).
option(OPENMS_WASM_THREADS "Build OpenMS for WebAssembly with threads (OpenMP, std::thread); needs the dependencies of the wasm32-emscripten-mt vcpkg triplet" OFF)

## 64-bit memory (Memory64, wasm64; CMake presets wasm64-emscripten-mt-release and
## wasm64-emscripten-mt-onnx-release): pointers and size_t are 64 bits, and the module can address
## more than the 4 GiB of wasm32 (browsers allow up to 16 GiB). Threaded only (the vcpkg triplet
## wasm64-emscripten-mt). -m64 (Emscripten deprecates -sMEMORY64) has to be in CMAKE_C_FLAGS and
## CMAKE_CXX_FLAGS when the toolchain is loaded (the presets set them), as Emscripten.cmake takes the
## pointer size from them.
option(OPENMS_WASM_MEMORY64 "Build OpenMS for 64-bit WebAssembly memory (Memory64); needs the dependencies of the wasm64-emscripten-mt vcpkg triplet" OFF)

## Not available in WebAssembly; openms_emscripten_check_options() rejects them.
set(_openms_emscripten_unsupported_options
  WITH_GUI              # Qt widgets
  WITH_THERMO_RAW       # needs the .NET runtime
  WITH_OPENTIMS         # Bruker .d directories are read with threads and sqlite's mmap
  WITH_CURL             # no sockets; a browser fetches through JavaScript
  ENABLE_UPDATE_CHECK   # needs WITH_CURL
  PYOPENMS
)
if(NOT OPENMS_WASM_THREADS)
  list(APPEND _openms_emscripten_unsupported_options
    MT_ENABLE_OPENMP    # Emscripten's libomp runs on threads: OPENMS_WASM_THREADS
  )
endif()
foreach(_option IN LISTS _openms_emscripten_unsupported_options)
  set(${_option} OFF CACHE BOOL "Not supported by the WebAssembly (Emscripten) build.")
endforeach()
unset(_option)

function(openms_emscripten_check_options)
  set(_enabled)
  foreach(_option IN LISTS _openms_emscripten_unsupported_options)
    if(${_option})
      list(APPEND _enabled ${_option})
    endif()
  endforeach()
  if(_enabled)
    list(JOIN _enabled ", " _enabled)
    message(FATAL_ERROR "Not supported by the WebAssembly (Emscripten) build: ${_enabled}. "
                        "Configure with these options OFF.")
  endif()
  ## The dependencies have to be built for the same threading mode as OpenMS.
  if(OPENMS_USE_VCPKG AND VCPKG_TARGET_TRIPLET)
    if(OPENMS_WASM_THREADS AND NOT VCPKG_TARGET_TRIPLET MATCHES "-mt$")
      message(FATAL_ERROR "OPENMS_WASM_THREADS needs dependencies built with threads: use the vcpkg "
                          "triplet wasm32-emscripten-mt (CMake preset wasm32-emscripten-mt-release), "
                          "not ${VCPKG_TARGET_TRIPLET}.")
    elseif(NOT OPENMS_WASM_THREADS AND VCPKG_TARGET_TRIPLET MATCHES "-mt$")
      message(FATAL_ERROR "The vcpkg triplet ${VCPKG_TARGET_TRIPLET} builds the dependencies with threads: "
                          "configure with OPENMS_WASM_THREADS=ON (CMake preset wasm32-emscripten-mt-release).")
    endif()
  endif()
  ## Memory64: with threads (the only wasm64 triplet, wasm64-emscripten-mt, and the presets), pointer size, and
  ## dependencies built for it
  if(OPENMS_WASM_MEMORY64 AND NOT OPENMS_WASM_THREADS)
    message(FATAL_ERROR "OPENMS_WASM_MEMORY64 needs OPENMS_WASM_THREADS (there is no single-threaded wasm64 build): use the "
                        "CMake preset wasm64-emscripten-mt-release.")
  endif()
  if(OPENMS_WASM_MEMORY64 AND NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(FATAL_ERROR "OPENMS_WASM_MEMORY64 needs -m64 in CMAKE_C_FLAGS and CMAKE_CXX_FLAGS "
                        "when the toolchain is loaded: use the CMake preset wasm64-emscripten-mt-release.")
  endif()
  if(OPENMS_USE_VCPKG AND VCPKG_TARGET_TRIPLET)
    if(OPENMS_WASM_MEMORY64 AND NOT VCPKG_TARGET_TRIPLET MATCHES "^wasm64")
      message(FATAL_ERROR "OPENMS_WASM_MEMORY64 needs dependencies built for 64-bit memory: use the vcpkg triplet "
                          "wasm64-emscripten-mt (CMake preset wasm64-emscripten-mt-release), not ${VCPKG_TARGET_TRIPLET}.")
    elseif(NOT OPENMS_WASM_MEMORY64 AND VCPKG_TARGET_TRIPLET MATCHES "^wasm64")
      message(FATAL_ERROR "The vcpkg triplet ${VCPKG_TARGET_TRIPLET} builds the dependencies for 64-bit memory: "
                          "configure with OPENMS_WASM_MEMORY64=ON (CMake preset wasm64-emscripten-mt-release).")
    endif()
  endif()
  ## The OpenMP default of the threaded build does not override a cached MT_ENABLE_OPENMP=OFF.
  if(OPENMS_WASM_THREADS AND NOT MT_ENABLE_OPENMP)
    message(FATAL_ERROR "OPENMS_WASM_THREADS needs OpenMP, but MT_ENABLE_OPENMP is OFF (usually kept in the cache "
                        "of a build directory configured without threads first): configure with "
                        "-DMT_ENABLE_OPENMP=ON, or use a fresh build directory.")
  endif()
endfunction()

#------------------------------------------------------------------------------
# Compiler and linker settings for every target of the build, the bundled
# third-party libraries and FetchContent dependencies included.

## OpenMS reports errors with C++ exceptions, so they have to be catchable: with
## Emscripten's default, a throw aborts the program. All objects linked into one
## module must agree on the exception handling mode; the wasm32-emscripten vcpkg
## triplet (vcpkg-overlays/triplets) builds the dependencies the same way.
## Native WebAssembly exception handling is supported by all current browsers and
## Node.js >= 17. The flag is also part of the usage requirements of the OpenMS
## targets (openms_add_compiler_flags), so code linking them inherits it.
set(OPENMS_EMSCRIPTEN_EH_FLAG "-fwasm-exceptions")
add_compile_options(${OPENMS_EMSCRIPTEN_EH_FLAG})
add_link_options(${OPENMS_EMSCRIPTEN_EH_FLAG})

## Threads: -pthread for every object (atomics, bulk memory, thread-local storage) and
## for the link (shared memory). Also part of the usage requirements of the OpenMS
## targets (openms_add_compiler_flags).
if(OPENMS_WASM_THREADS)
  add_compile_options(-pthread)
  add_link_options(-pthread)
  ## OpenMP. The runtime is a system library of Emscripten, linked by -fopenmp on the
  ## link line (the compile flag comes with OpenMP::OpenMP_CXX). Preset for FindOpenMP:
  ## it would otherwise record the system libraries of its test program (debug
  ## variants of libc, libc++ and malloc) as link libraries of OpenMP::OpenMP_CXX.
  add_link_options(-fopenmp)
  set(MT_ENABLE_OPENMP ON CACHE BOOL "Enable OpenMP support")
  set(MT_REQUIRE_OPENMP ON CACHE BOOL "With MT_ENABLE_OPENMP, fail the configuration if no OpenMP runtime is found, instead of falling back to -fopenmp-simd (no parallel regions). The CI presets set it.")
  foreach(_lang C CXX)
    set(OpenMP_${_lang}_FLAGS "-fopenmp" CACHE STRING "Compiler flag of OpenMP (Emscripten)")
    set(OpenMP_${_lang}_LIB_NAMES "" CACHE STRING "Libraries of OpenMP (Emscripten: none, -fopenmp links the runtime)")
  endforeach()
  unset(_lang)
  ## No nested parallelism: every thread of a team is a Web Worker, and in a browser the
  ## workers have to be started before a team needs them (src/openms_wasm sizes that pool).
  set(MT_ENABLE_NESTED_OPENMP OFF CACHE BOOL "Enable nested parallelism.")
endif()

## mimalloc instead of Emscripten's default allocator (dlmalloc). With threads, dlmalloc
## serializes all of them on one lock: a database search (SimpleSearchEngine) ran 1.8x slower
## on 2 threads and 3.4x slower on 4 threads than on 1; with mimalloc it is 1.4x and 1.7x
## faster. Single-threaded it is about as fast or faster: search 1549 instead of 1601 ms,
## loading a 14 MB mzML file 184 instead of 238 ms. The module grows by about 0.3 MB.
## OpenMS' build of it (OPENMS_WASM_MALLOC, the default) or Emscripten's: see
## cmake/emscripten_malloc.cmake.
include(${CMAKE_CURRENT_LIST_DIR}/emscripten_malloc.cmake)

## The JavaScript bindings (src/openms_wasm): the module web applications load.
option(OPENMS_WASM_BINDINGS "Build the JavaScript bindings of the OpenMS library (src/openms_wasm)" ON)

## Where the library looks for its shared data first (share/OpenMS): the JavaScript module
## preloads it into this directory of its virtual file system. (Programs run with Node.js,
## like the class tests, find it in the source tree instead.)
set(OPENMS_INSTALL_DATA_PATH "/share/OpenMS" CACHE PATH "Compiled-in location of the installed shared data (share/OpenMS) of OpenMS. Empty: <CMAKE_INSTALL_PREFIX>/<INSTALL_SHARE_DIR>.")

## Fixed-width SIMD (128 bit) is supported by all current browsers and Node.js >= 16.4;
## OpenMS' SSE code paths (via SIMDe) and the compiler's auto-vectorizer use it.
option(OPENMS_WASM_SIMD "Compile OpenMS with WebAssembly SIMD (-msimd128)" ON)

## Memory and stack of every executable (the JavaScript module, and the TOPP tools
## and class tests where they are built). Emscripten's default stack (64 KiB) is
## far smaller than the 8 MiB of a native main thread, which OpenMS (recursive
## algorithms, large objects on the stack, Xerces) relies on. The heap starts at
## 64 MiB (room for the stack and the shared data the module embeds) and grows on
## demand up to the 4 GiB a wasm32 module can address, or with 64-bit memory up to
## 16 GiB, the most browsers (and Node.js) give a Memory64 module.
set(OPENMS_WASM_STACK_SIZE "8MB" CACHE STRING "Stack size of OpenMS WebAssembly executables (Emscripten STACK_SIZE)")
if(OPENMS_WASM_MEMORY64)
  add_link_options("-m64")
  set(_openms_wasm_maximum_memory "16GB")
else()
  set(_openms_wasm_maximum_memory "4GB")
endif()
add_link_options(
  "-sALLOW_MEMORY_GROWTH=1"
  "-sINITIAL_MEMORY=64MB"
  "-sMAXIMUM_MEMORY=${_openms_wasm_maximum_memory}"
  "-sSTACK_SIZE=${OPENMS_WASM_STACK_SIZE}"
)

## Link options of the programs OpenMS runs with Node.js (class tests, TOPP tools):
## they read and write the files of the host directly (NODERAWFS) and return their
## exit code (EXIT_RUNTIME). EXECUTABLE starts the JavaScript file with a
## '#!/usr/bin/env node' line and makes it executable.
set(OPENMS_EMSCRIPTEN_NODE_PROGRAM_LINK_OPTIONS
  "-sNODERAWFS=1"
  "-sEXIT_RUNTIME=1"
  "-sENVIRONMENT=node"
  "-sEXECUTABLE=1"
)

## Makes executable @p target such a Node.js program (openms_add_executable_compiler_flags
## calls it for the executables of OpenMS; targets that do not use that function call it
## themselves). It also drops the '.js' suffix Emscripten's toolchain file gives executables
## (per target: project() calls of subdirectories reset CMAKE_EXECUTABLE_SUFFIX): <bin>/Foo
## then runs like a native program, so the tests run it the same way on every platform.
function(openms_emscripten_node_program target)
  target_link_options(${target} PRIVATE ${OPENMS_EMSCRIPTEN_NODE_PROGRAM_LINK_OPTIONS})
  openms_emscripten_link_malloc(${target})
  set_target_properties(${target} PROPERTIES SUFFIX "")
endfunction()

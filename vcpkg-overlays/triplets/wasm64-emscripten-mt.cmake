# WebAssembly with 64-bit memory (Memory64, wasm64) for the multi-threaded Emscripten build of the
# OpenMS library (CMake presets wasm64-emscripten-mt-release and wasm64-emscripten-mt-onnx-release,
# see src/openms_wasm/README.md).
#
# wasm32-emscripten-mt.cmake plus -m64 (Memory64): pointers and size_t are 64 bits, and the module
# can address more than 4 GiB. Keep the files in sync (vcpkg hashes the triplet file itself, so
# this one does not include the other). vcpkg has no wasm64 architecture: the architecture stays
# wasm32 (a label), the flags and the toolchain (toolchains/wasm64-emscripten.cmake, which sets
# the pointer size) make the difference.
#
# wasm32-emscripten-mt.cmake is wasm32-emscripten.cmake plus -pthread: a module with threads
# (shared memory) can only be linked from objects compiled for it, so the threaded build needs
# its own build of every dependency.
#
# It replaces vcpkg's community triplet of the same name. The difference: every
# dependency is compiled with native WebAssembly exception handling
# (-fwasm-exceptions). OpenMS reports errors through C++ exceptions, so catching
# has to work, and Emscripten requires one exception handling mode for all objects
# of a program; the community triplet builds with Emscripten's default, which
# aborts on any throw. Native Wasm exceptions are supported by every current
# browser and Node.js (>= 17) and are faster and smaller than the JavaScript-based
# emulation (-fexceptions).
#
# Emscripten is located like the community triplet does it: EMSCRIPTEN_ROOT, emcc
# on PATH (source emsdk_env.sh) or EMSDK.

set(VCPKG_ENV_PASSTHROUGH_UNTRACKED EMSCRIPTEN_ROOT EMSDK PATH)

if(NOT DEFINED ENV{EMSCRIPTEN_ROOT})
  find_path(EMSCRIPTEN_ROOT "emcc")
else()
  set(EMSCRIPTEN_ROOT "$ENV{EMSCRIPTEN_ROOT}")
endif()

if(NOT EMSCRIPTEN_ROOT)
  if(NOT DEFINED ENV{EMSDK})
    message(FATAL_ERROR "The emcc compiler was not found on PATH; source emsdk_env.sh or set EMSDK.")
  endif()
  set(EMSCRIPTEN_ROOT "$ENV{EMSDK}/upstream/emscripten")
endif()

if(NOT EXISTS "${EMSCRIPTEN_ROOT}/cmake/Modules/Platform/Emscripten.cmake")
  message(FATAL_ERROR "Emscripten.cmake toolchain file not found below ${EMSCRIPTEN_ROOT}")
endif()

set(VCPKG_TARGET_ARCHITECTURE wasm32)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE static)
## OpenMS builds only Release for WebAssembly; skip the debug half of every port.
set(VCPKG_BUILD_TYPE release)
set(VCPKG_CMAKE_SYSTEM_NAME Emscripten)
## Emscripten.cmake plus the flags below: vcpkg applies VCPKG_*_FLAGS only through its
## own platform toolchains, not through a chainloaded one (see the toolchain file).
get_filename_component(VCPKG_CHAINLOAD_TOOLCHAIN_FILE
  "${CMAKE_CURRENT_LIST_DIR}/../toolchains/wasm64-emscripten.cmake" ABSOLUTE)

## C as well: setjmp/longjmp support follows the exception handling mode, and the
## two must agree between the objects linked into one module. -pthread: atomics and
## bulk memory, thread-local storage. -m64: 64-bit memory (wasm64; Emscripten deprecates -sMEMORY64).
set(VCPKG_C_FLAGS "-fwasm-exceptions -pthread -m64")
set(VCPKG_CXX_FLAGS "-fwasm-exceptions -pthread -m64")
set(VCPKG_LINKER_FLAGS "-fwasm-exceptions -pthread -m64")

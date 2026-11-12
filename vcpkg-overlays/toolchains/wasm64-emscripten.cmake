# Chainloaded toolchain of the wasm64-emscripten-mt overlay triplet
# (vcpkg-overlays/triplets/wasm64-emscripten-mt.cmake): wasm32-emscripten.cmake for 64-bit
# memory (Memory64, -m64 in the triplet's flags).
#
# A copy of wasm32-emscripten.cmake plus the pointer size; keep the two in sync. (vcpkg
# hashes the chainloaded toolchain file itself into the ABI of every port, so this one does
# not include the other, and the wasm32 triplets keep their binary caches.)
#
# Emscripten.cmake takes the pointer size (CMAKE_SIZEOF_VOID_P, CMAKE_LIBRARY_ARCHITECTURE)
# from CMAKE_C_FLAGS. CMake 4 loads it a second time, as platform file in project(), before
# it initializes CMAKE_C_FLAGS from CMAKE_C_FLAGS_INIT, which then holds the triplet's
# flags: it would take the module for wasm32. So the triplet's flags go into the cache.

if(NOT _OPENMS_WASM64_EMSCRIPTEN_TOOLCHAIN)
  set(_OPENMS_WASM64_EMSCRIPTEN_TOOLCHAIN 1)

  # Emscripten is located like the triplet does it (the environment is passed through).
  if(DEFINED ENV{EMSCRIPTEN_ROOT})
    set(_openms_emscripten_root "$ENV{EMSCRIPTEN_ROOT}")
  else()
    find_path(_openms_emscripten_root "emcc")
    if(NOT _openms_emscripten_root AND DEFINED ENV{EMSDK})
      set(_openms_emscripten_root "$ENV{EMSDK}/upstream/emscripten")
    endif()
  endif()
  if(NOT EXISTS "${_openms_emscripten_root}/cmake/Modules/Platform/Emscripten.cmake")
    message(FATAL_ERROR "Emscripten.cmake toolchain file not found; source emsdk_env.sh or set EMSDK.")
  endif()
  if(NOT " ${VCPKG_C_FLAGS} " MATCHES " -m64 ")
    message(FATAL_ERROR "wasm64-emscripten.cmake is the toolchain of a triplet with -m64 in VCPKG_C_FLAGS.")
  endif()
  # The triplet's flags into the cache (see above); CMAKE_*_FLAGS_INIT below are then not used.
  if(NOT DEFINED CACHE{CMAKE_C_FLAGS})
    set(CMAKE_C_FLAGS "${VCPKG_C_FLAGS}" CACHE STRING "Flags used by the C compiler (wasm64 triplet)")
  endif()
  if(NOT DEFINED CACHE{CMAKE_CXX_FLAGS})
    set(CMAKE_CXX_FLAGS "${VCPKG_CXX_FLAGS}" CACHE STRING "Flags used by the C++ compiler (wasm64 triplet)")
  endif()
  include("${_openms_emscripten_root}/cmake/Modules/Platform/Emscripten.cmake")

  # Also where a port sets CMAKE_C_FLAGS itself (its flags get the triplet's appended below).
  # (CMake derives CMAKE_SIZEOF_VOID_P from the compilers' data pointer size again: set all three.)
  set(CMAKE_LIBRARY_ARCHITECTURE "wasm64-emscripten")
  set(CMAKE_SIZEOF_VOID_P 8)
  set(CMAKE_C_SIZEOF_DATA_PTR 8)
  set(CMAKE_CXX_SIZEOF_DATA_PTR 8)

  # As in vcpkg's platform toolchains: try_compile() projects see the same flags.
  list(APPEND CMAKE_TRY_COMPILE_PLATFORM_VARIABLES
    VCPKG_C_FLAGS VCPKG_CXX_FLAGS
    VCPKG_C_FLAGS_DEBUG VCPKG_CXX_FLAGS_DEBUG
    VCPKG_C_FLAGS_RELEASE VCPKG_CXX_FLAGS_RELEASE
    VCPKG_LINKER_FLAGS VCPKG_LINKER_FLAGS_RELEASE VCPKG_LINKER_FLAGS_DEBUG
  )

  string(APPEND CMAKE_C_FLAGS_INIT " ${VCPKG_C_FLAGS} ")
  string(APPEND CMAKE_CXX_FLAGS_INIT " ${VCPKG_CXX_FLAGS} ")
  string(APPEND CMAKE_C_FLAGS_DEBUG_INIT " ${VCPKG_C_FLAGS_DEBUG} ")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG_INIT " ${VCPKG_CXX_FLAGS_DEBUG} ")
  string(APPEND CMAKE_C_FLAGS_RELEASE_INIT " ${VCPKG_C_FLAGS_RELEASE} ")
  string(APPEND CMAKE_CXX_FLAGS_RELEASE_INIT " ${VCPKG_CXX_FLAGS_RELEASE} ")

  string(APPEND CMAKE_MODULE_LINKER_FLAGS_INIT " ${VCPKG_LINKER_FLAGS} ")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS_INIT " ${VCPKG_LINKER_FLAGS} ")
  string(APPEND CMAKE_EXE_LINKER_FLAGS_INIT " ${VCPKG_LINKER_FLAGS} ")
  string(APPEND CMAKE_MODULE_LINKER_FLAGS_DEBUG_INIT " ${VCPKG_LINKER_FLAGS_DEBUG} ")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS_DEBUG_INIT " ${VCPKG_LINKER_FLAGS_DEBUG} ")
  string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG_INIT " ${VCPKG_LINKER_FLAGS_DEBUG} ")
  string(APPEND CMAKE_MODULE_LINKER_FLAGS_RELEASE_INIT " ${VCPKG_LINKER_FLAGS_RELEASE} ")
  string(APPEND CMAKE_SHARED_LINKER_FLAGS_RELEASE_INIT " ${VCPKG_LINKER_FLAGS_RELEASE} ")
  string(APPEND CMAKE_EXE_LINKER_FLAGS_RELEASE_INIT " ${VCPKG_LINKER_FLAGS_RELEASE} ")
endif()

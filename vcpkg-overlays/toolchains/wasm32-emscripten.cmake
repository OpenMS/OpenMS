# Chainloaded toolchain of the wasm32-emscripten overlay triplet
# (vcpkg-overlays/triplets/wasm32-emscripten.cmake): Emscripten's toolchain file plus
# the compiler and linker flags of the triplet.
#
# vcpkg applies VCPKG_C_FLAGS, VCPKG_CXX_FLAGS and VCPKG_LINKER_FLAGS only through its
# own platform toolchains (scripts/toolchains/linux.cmake, ...); a triplet that
# chainloads Emscripten.cmake directly builds every CMake port without them. Then the
# dependencies lack -fwasm-exceptions: their catch blocks are compiled out, and an
# exception they throw and catch internally (Xerces-C does so for every recoverable
# parse error) escapes into the caller instead.

if(NOT _OPENMS_WASM32_EMSCRIPTEN_TOOLCHAIN)
  set(_OPENMS_WASM32_EMSCRIPTEN_TOOLCHAIN 1)

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
  include("${_openms_emscripten_root}/cmake/Modules/Platform/Emscripten.cmake")

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

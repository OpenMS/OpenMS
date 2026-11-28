# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# The allocator of OpenMS' WebAssembly programs (included by cmake/emscripten.cmake).
#
# OPENMS_WASM_MALLOC:
#  mimalloc          Emscripten's mimalloc (-sMALLOC=mimalloc).
#  openms-mimalloc   the same sources of mimalloc (those of the Emscripten SDK), built here for the
#                    programs that link the target openms_wasm_malloc (the JavaScript module, and the
#                    Node.js programs of openms_emscripten_node_program()), with -sMALLOC=none:
#                    - MI_FREE_USE_PAGEMAP: free() finds the page of a block in mimalloc's page map, so
#                      that the page metadata does not have to be aligned. Otherwise mimalloc aligns its
#                      arenas to 64 MiB in wasm64 (16 MiB in wasm32) and takes every block larger than
#                      that from the system allocator (below) instead of from an arena.
#                    - no exponential growth of the arenas (mimalloc doubles the size of the arena it
#                      reserves every 8 arenas: GBs per arena in a long wasm64 search).
#                    - dlmalloc as the system allocator below mimalloc instead of emmalloc
#                      (src/openms_wasm/malloc/dlmalloc_system.c), which in wasm64 cannot allocate a
#                      block of 4 GiB or more, corrupts free regions of 4 GiB or more
#                      (https://github.com/emscripten-core/emscripten/issues/28028), and reuses freed
#                      huge blocks poorly.
#                    The default of every build, wasm32 and wasm64 (an existing build directory keeps
#                    the value in its cache). ProSE, the memory of the module after a search
#                    (getMemoryInfo().heapSize), as fast with both; wasm64: Node.js, wasm32: Chromium 141
#                    (searchFiles(), file by file, as the ProSE app), 16 threads:
#                    | search                                  | wasm64: Emscripten's, this | wasm32: Emscripten's, this |
#                    | Erwinia, tryptic                        | 1.9 GB, 0.95 GB            | 1.0 GB, 0.67 GB            |
#                    | human, tryptic                          | 2.5 GB, 1.35 GB            | 1.5 GB, 0.94 GB            |
#                    | human, semi-tryptic, chunk_size 2000    | 12.2 GB, 1.75 GB           | 1.9-3.9 GB, 1.5-1.7 GB     |
#                    | human, semi-tryptic, chunk_size 8000    | 10.5 GB, 4.7 GB            | (does not fit)             |
#                    | 3 x Erwinia, tryptic, file by file      |                            | 1.5 GB, 0.8 GB             |
#                    wasm32 single-threaded: Erwinia 0.59 and 0.39 GB, human 0.99 and 0.79 GB.
#------------------------------------------------------------------------------

set(OPENMS_WASM_MALLOC "openms-mimalloc" CACHE STRING
    "Allocator of the OpenMS WebAssembly programs: openms-mimalloc (mimalloc built by OpenMS, see cmake/emscripten_malloc.cmake; the default) or mimalloc (Emscripten's)")
set_property(CACHE OPENMS_WASM_MALLOC PROPERTY STRINGS mimalloc openms-mimalloc)
if(NOT OPENMS_WASM_MALLOC MATCHES "^(mimalloc|openms-mimalloc)$")
  message(FATAL_ERROR "OPENMS_WASM_MALLOC must be mimalloc or openms-mimalloc (is '${OPENMS_WASM_MALLOC}').")
endif()

## The configuration checks that the sources are as this file expects them (mimalloc, dlmalloc, mimalloc's calls of
## emmalloc_memalign()/emmalloc_free(), the exports emscripten_builtin_* of dlmalloc.c, and the growth of the arenas in
## mimalloc's arena.c): if a new Emscripten changed them, it fails with a hint.

## Every program: Emscripten's mimalloc, or none where the program links openms_wasm_malloc (openms_emscripten_link_malloc()
## sets the target property OPENMS_WASM_OWN_MALLOC): one -sMALLOC per link, whatever the order of the options.
add_link_options("-sMALLOC=$<IF:$<BOOL:$<TARGET_PROPERTY:OPENMS_WASM_OWN_MALLOC>>,none,mimalloc>")

if(OPENMS_WASM_MALLOC STREQUAL "openms-mimalloc")
  enable_language(C)
  set(_mi "${EMSCRIPTEN_ROOT_PATH}/system/lib/mimalloc")
  if(NOT EXISTS "${_mi}/src/prim/emscripten/prim.c" OR NOT EXISTS "${EMSCRIPTEN_ROOT_PATH}/system/lib/dlmalloc.c")
    message(FATAL_ERROR "OPENMS_WASM_MALLOC=openms-mimalloc: the sources of mimalloc or dlmalloc are not where Emscripten "
                        "${EMSCRIPTEN_VERSION} had them (${_mi}); use OPENMS_WASM_MALLOC=mimalloc, or update cmake/emscripten_malloc.cmake.")
  endif()
  file(STRINGS "${_mi}/src/prim/emscripten/prim.c" _calls REGEX "emmalloc_(memalign|free)\\(")
  if(NOT _calls)
    message(FATAL_ERROR "OPENMS_WASM_MALLOC=openms-mimalloc: mimalloc of Emscripten ${EMSCRIPTEN_VERSION} no longer takes its memory "
                        "from emmalloc_memalign() and emmalloc_free(); update src/openms_wasm/malloc/dlmalloc_system.c.")
  endif()
  ## dlmalloc.c defines these emscripten_builtin_* (whatever USE_DL_PREFIX says), which dlmalloc_system.c renames: another
  ## one would collide with mimalloc's, or replace it
  file(READ "${EMSCRIPTEN_ROOT_PATH}/system/lib/dlmalloc.c" _dlmalloc)
  string(REGEX MATCHALL "emscripten_builtin_[a-z_]+" _builtins "${_dlmalloc}")
  list(REMOVE_DUPLICATES _builtins)
  list(SORT _builtins)
  if(NOT _builtins STREQUAL "emscripten_builtin_calloc;emscripten_builtin_free;emscripten_builtin_malloc;emscripten_builtin_memalign;emscripten_builtin_realloc")
    message(FATAL_ERROR "OPENMS_WASM_MALLOC=openms-mimalloc: dlmalloc.c of Emscripten ${EMSCRIPTEN_VERSION} defines ${_builtins}; "
                        "update the renames of src/openms_wasm/malloc/dlmalloc_system.c and this check (cmake/emscripten_malloc.cmake).")
  endif()

  ## The sources Emscripten builds libmimalloc from (tools/system_libs.py), but emmalloc.
  file(GLOB _mi_sources "${_mi}/src/*.c")
  list(FILTER _mi_sources EXCLUDE REGEX "/(alloc-override|free|page-queue|static|arena)\\.c$") # (included by others; arena.c below)
  list(APPEND _mi_sources "${_mi}/src/prim/prim.c" "${_mi}/src/prim/prim-tls.c" "${EMSCRIPTEN_ROOT_PATH}/system/lib/libc/sbrk.c")

  ## arena.c without the exponential growth of the arenas (a copy; OPENMS_MI_ARENA_GROWTH)
  file(READ "${_mi}/src/arena.c" _arena)
  set(_growth "  if (arena_count >= 1 && arena_count <= 128) {\n    // scale up the arena sizes exponentially every 8 entries")
  string(FIND "${_arena}" "${_growth}" _at)
  ## (A failure, not a warning: with the growth, a wasm64 search took twice the memory; a new Emscripten's arena.c is to
  ## be looked at anyway.)
  if(_at EQUAL -1)
    message(FATAL_ERROR "OPENMS_WASM_MALLOC=openms-mimalloc: mimalloc's arena.c of Emscripten ${EMSCRIPTEN_VERSION} has changed: "
                        "update the patch of its exponential growth of arenas in cmake/emscripten_malloc.cmake, or use "
                        "OPENMS_WASM_MALLOC=mimalloc.")
  endif()
  string(REPLACE "${_growth}"
    "  if (OPENMS_MI_ARENA_GROWTH && arena_count >= 1 && arena_count <= 128) { // OpenMS: not by default (cmake/emscripten_malloc.cmake)\n    // scale up the arena sizes exponentially every 8 entries"
    _arena "${_arena}")
  file(WRITE "${CMAKE_BINARY_DIR}/openms_wasm_malloc/arena.c.new" "${_arena}")
  configure_file("${CMAKE_BINARY_DIR}/openms_wasm_malloc/arena.c.new" "${CMAKE_BINARY_DIR}/openms_wasm_malloc/arena.c" COPYONLY)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_mi}/src/arena.c")

  add_library(openms_wasm_malloc OBJECT
    ${_mi_sources}
    "${CMAKE_BINARY_DIR}/openms_wasm_malloc/arena.c"
    "${PROJECT_SOURCE_DIR}/src/openms_wasm/malloc/dlmalloc_system.c")
  target_include_directories(openms_wasm_malloc PRIVATE "${_mi}/include" "${_mi}/src" "${EMSCRIPTEN_ROOT_PATH}/system/lib")
  ## Emscripten's flags of libmimalloc, and OpenMS' changes
  target_compile_definitions(openms_wasm_malloc PRIVATE
    MI_ENABLE_LARGE_PAGES=0 "MI_ARENA_SLICE_SHIFT=(12+MI_SIZE_SHIFT)" MI_MAX_ALIGN_SIZE=8 MI_DEFAULT_ARENA_RESERVE=65536
    MI_MALLOC_OVERRIDE MI_DEBUG=0 NDEBUG MI_LIBC_MUSL MI_USE_BUILTIN_THREAD_POINTER
    MI_FREE_USE_PAGEMAP=1 OPENMS_MI_ARENA_GROWTH=0)
  target_compile_options(openms_wasm_malloc PRIVATE -fno-builtin -w)
endif()

## Links the allocator of OPENMS_WASM_MALLOC into the program @p target (the JavaScript module and the Node.js programs):
## with -sMALLOC=none (OPENMS_WASM_OWN_MALLOC, above).
function(openms_emscripten_link_malloc target)
  if(TARGET openms_wasm_malloc)
    target_link_libraries(${target} PRIVATE openms_wasm_malloc)
    set_target_properties(${target} PROPERTIES OPENMS_WASM_OWN_MALLOC TRUE)
  endif()
endfunction()

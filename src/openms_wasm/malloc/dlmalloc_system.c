// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// The system allocator under mimalloc in OpenMS' build of it (OPENMS_WASM_MALLOC=openms-mimalloc,
// cmake/emscripten_malloc.cmake): Emscripten's dlmalloc instead of emmalloc.
//
// mimalloc takes its arenas and its huge blocks (beyond its arenas) from the functions
// emmalloc_memalign() and emmalloc_free() (src/prim/emscripten/prim.c of Emscripten's mimalloc).
// Emscripten links emmalloc there, which in a wasm64 module (Emscripten 6.0.12) cannot allocate a
// block of 4 GiB or more (MAX_ALLOC_SIZE) and corrupts free regions of 4 GiB or more (their size is
// masked with a 32-bit ~FREE_REGION_FLAG, their free list computed from __builtin_clz of a size_t:
// https://github.com/emscripten-core/emscripten/issues/28028), and which reuses freed huge blocks
// poorly: a search that frees and allocates blocks of hundreds of MB again and again (ProSE
// with database:chunk_size) grew the memory of the wasm64 module several times beyond what it used.
// dlmalloc is 64-bit clean and coalesces freed blocks; its lock costs nothing here, as mimalloc
// calls it rarely.
//
// dlmalloc.c of the Emscripten SDK is compiled here as it is. Its exports emscripten_builtin_*
// would collide with mimalloc's and are renamed; its weak aliases of malloc, free, ... lose to
// mimalloc's definitions.
//
// mimalloc asks for blocks aligned to its slices (32 KiB in wasm64), mostly large ones. dlmemalign()
// serves such a request from a free block only if that is larger than the request by the alignment:
// a block of n bytes, once freed, was not given to the next request of n bytes with that alignment.
// A module that allocates and frees large blocks again and again grew by each, e.g. by the size of a
// file of the file system, whose memory allocateFile() probes before WasmFS allocates it (a file of
// 4 GiB after one of 8 GiB did not fit into 16 GiB, and WasmFS aborted). So the block is allocated
// with room for the alignment by dlmalloc(), which reuses a freed block of the same size, and its
// start is kept in front of the aligned address (at most the alignment and 8 bytes more per block;
// mimalloc allocates a few hundred blocks here during a large search).

#define emscripten_builtin_malloc openms_system_builtin_malloc
#define emscripten_builtin_realloc openms_system_builtin_realloc
#define emscripten_builtin_calloc openms_system_builtin_calloc
#define emscripten_builtin_free openms_system_builtin_free
#define emscripten_builtin_memalign openms_system_builtin_memalign

#include "dlmalloc.c" // system/lib of the Emscripten SDK

void* emmalloc_memalign(size_t alignment, size_t size) // (alignment: a power of 2)
{
  if (alignment < sizeof(void*)) alignment = sizeof(void*);
  if (size > SIZE_MAX - alignment - sizeof(void*)) return NULL;
  char* start = (char*)dlmalloc(size + alignment + sizeof(void*));
  if (!start) return NULL;
  uintptr_t block = ((uintptr_t)start + sizeof(void*) + alignment - 1) & ~(uintptr_t)(alignment - 1);
  ((void**)block)[-1] = start;
  return (void*)block;
}

void emmalloc_free(void* block)
{
  if (block) dlfree(((void**)block)[-1]);
}

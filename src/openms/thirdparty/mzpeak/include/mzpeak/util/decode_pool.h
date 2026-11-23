/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <arrow/memory_pool.h>

#include <atomic>
#include <cstdint>
#include <mutex>
#include <string>

namespace MzPeak::Util {

/******************************************************************************/
/// The Arrow memory pool Parquet decoding allocates from: the decoded arrays,
/// Parquet's page and decompression buffers, and the column-chunk reads under
/// them.
///
/// WHY NOT MALLOC.  A decoded row group is freed by whichever reader evicts
/// it from the shared RowGroupCache, which is rarely the thread that decoded
/// it.  glibc returns a freed chunk to the arena of the thread that ALLOCATED
/// it, so every decoding thread's arena ends up holding about one freed group
/// (~26 MB) that only that thread would ever reuse: 3.56 GB in 137 arenas at
/// 128 threads, and RSS that grows with the thread count.  Arrow's bundled
/// jemalloc and mimalloc are no way out either: they have crashed on a buffer
/// freed after the thread that allocated it had exited (see Parquet::Impl),
/// and a consumer's OpenMP or std::async threads come and go.
///
/// WHAT IT DOES.  Buffers of kMinBlock and more are whole blocks mapped from
/// the OS, their size rounded up to a size class: four per power of two, but
/// never finer than the 64 KB granule.  So a block is under 25% over its
/// request from 256 KB up, under 50% for 128-256 KB and under 100% for
/// 64-128 KB (65,537 bytes take 128 KB).  A page takes memory once it is
/// touched (on Linux, once it is written), so on macOS a fresh block's slack
/// costs address space, not memory.  On Linux the same holds for small pages,
/// but slack that shares a huge page with written bytes is resident with it
/// (see below).  On Windows the whole block is committed when it is mapped,
/// so its slack, like a cached block, counts against the commit limit (RAM
/// plus page file) from the start, though its pages take RAM only once
/// touched.  A freed block goes onto its class's free list, shared by all
/// threads, and the next request of that class on any thread takes the most
/// recently freed one.  The block keeps its pages until it is unmapped, so a
/// smaller request of its class that reuses it carries the pages its earlier
/// users wrote: in steady state the pool's memory is about bytes_mapped(), not
/// the bytes requested.  Which thread frees a block, and whether the thread
/// that mapped it still exists, does not matter.  Smaller buffers, and any
/// alignment above 64, go to Arrow's system allocator as before.
///
/// On Linux the blocks are backed by transparent huge pages: smaller ones are
/// carved from 64 MB regions, and blocks of 2 MB and more are mapped on their
/// own, both 2 MB-aligned and advised for huge pages.  A fresh 4 KB page costs
/// ~10 us to fault in with a hundred threads running, and 4 KB faults were
/// what separated this pool from Arrow's jemalloc on speed: the ~4 MB column
/// chunks Arrow's eight I/O threads read for every reader took 8x longer into
/// fresh 4 KB pages than into warm memory.  Measured at 128 threads: minor
/// faults 1.0 -> 0.5 million, wall 3.8 -> 3.3 s, same RSS.  The price is RSS
/// in 2 MB steps: a block's untouched class slack is resident once its huge
/// page is.
///
/// A block's size follows from the size Arrow passes to Free() and
/// Reallocate(), which Arrow's MemoryPool contract makes the size last
/// allocated (jemalloc's sized deallocation depends on the same); the pool
/// keeps no per-block record.
///
/// THE BOUND.  Freed blocks are kept only while
///
///   cached <= max(retain_floor, peak / 4, peak - live)
///
/// where live is the bytes of allocated blocks and peak its high-water mark
/// over the current and the previous second; the oldest-freed blocks beyond
/// that are unmapped.  So the pool's footprint (live + cached) stays within
/// its own recent peak plus a quarter: memory an eviction burst releases is
/// kept for the decodes the burst made room for, the quarter absorbs the
/// churn of blocks whose size class does not come back at once, and nothing
/// goes to an arena that strands it.  A peak that has passed is released on
/// the next allocation or free two seconds on.
///
/// Chosen over a fixed-cap free list and over per-row-group slabs by replaying
/// a recorded allocation trace (128 threads, 347 groups); the numbers are in
/// the commit message.
class DecodePool final : public arrow::MemoryPool {
public:
  /// Smallest buffer served from a block; smaller ones go to the system pool.
  static constexpr int64_t kMinBlock = int64_t{64} << 10;

  /// Freed bytes that may always be kept, however small the peak: a few
  /// decodes' worth of page and decompression buffers.  ponytail: a
  /// constant; above it the bound scales with the peak.
  static constexpr int64_t kDefaultRetainFloor = int64_t{64} << 20;

  /// The huge page size; blocks this size and up get their own mapping.
  static constexpr int64_t kHugePage = int64_t{2} << 20;

  explicit DecodePool(int64_t retain_floor = kDefaultRetainFloor);

  /// Unmaps the cached blocks and the uncarved rest of the current region.
  /// Blocks still allocated stay mapped.
  ~DecodePool() override;

  DecodePool(const DecodePool&) = delete;
  DecodePool& operator=(const DecodePool&) = delete;

  // The default-alignment overloads, which the overrides would hide.
  using arrow::MemoryPool::Allocate;
  using arrow::MemoryPool::Free;
  using arrow::MemoryPool::Reallocate;

  arrow::Status Allocate(int64_t size, int64_t alignment, uint8_t** out) override;
  arrow::Status Reallocate(int64_t old_size,
                           int64_t new_size,
                           int64_t alignment,
                           uint8_t** ptr) override;
  void Free(uint8_t* buffer, int64_t size, int64_t alignment) override;

  /// Unmaps every cached block.
  void ReleaseUnused() override;

  /// Arrow's accounting, over everything allocated through this pool (blocks
  /// and system-pool buffers alike), in requested bytes.
  int64_t bytes_allocated() const override;
  int64_t max_memory() const override;
  int64_t total_bytes_allocated() const override;
  int64_t num_allocations() const override;
  std::string backend_name() const override;

  /// Bytes of blocks, allocated and cached.  The pool's share of RSS is about
  /// this: less by what was never written, more by huge-page rounding.  On
  /// Windows it is also the blocks' commit charge.
  int64_t bytes_mapped() const;

  /// Bytes of freed blocks kept for reuse.
  int64_t bytes_cached() const;

  /// The size of the block a request of @p size is served from.
  static int64_t block_size(int64_t size);

private:
  struct Node; // a cached block's own first bytes

  static constexpr int kClasses = 256;
  static constexpr int64_t kMaxAlignment = 64;
  static constexpr int64_t kRegion = int64_t{64} << 20;
  // ponytail: bigger requests go to the system pool; nothing decodes 256 TB.
  static constexpr int64_t kMaxBlock = int64_t{1} << 48;

  static bool block_backed_(int64_t size, int64_t alignment)
  {
    return size >= kMinBlock && size <= kMaxBlock && alignment <= kMaxAlignment;
  }

  /// A block of @p block bytes, cached or freshly mapped; null when out of memory.
  uint8_t* take_(int64_t block);
  /// Caller holds mutex_.  A fresh block from the current region (Linux).
  uint8_t* carve_(int64_t block);
  /// Caches @p p, then unmaps the oldest cached blocks beyond the bound.
  void give_(uint8_t* p, int64_t block);
  /// Caller holds mutex_.  Rolls the one-second epochs of the live peak.
  void roll_epoch_();
  /// Caller holds mutex_.  Unlinks @p n from both lists.
  void unlink_(Node* n);
  /// Caller holds mutex_.  Unlinks the oldest cached blocks beyond the bound
  /// and returns them, chained through Node::next, for unmap_().
  Node* trim_();
  /// Unmaps a chain from trim_(); called without the lock.
  static void unmap_(Node* doomed);

  void did_allocate_(int64_t bytes);
  void did_free_(int64_t bytes);

  arrow::MemoryPool* const system_;
  const int64_t retain_floor_;

  mutable std::mutex mutex_;
  Node* lists_[kClasses] = {}; ///< per class, most recently freed first
  Node* oldest_ = nullptr;     ///< across classes, for unmapping
  Node* newest_ = nullptr;
  int64_t live_ = 0;   ///< mapped bytes of allocated blocks
  int64_t cached_ = 0; ///< mapped bytes of cached blocks
  int64_t epoch_ = 0;
  int64_t peak_ = 0;      ///< largest live_ this epoch
  int64_t last_peak_ = 0; ///< ... and the one before
  uint8_t* region_next_ = nullptr; ///< carve_()'s current region
  uint8_t* region_end_ = nullptr;

  std::atomic<int64_t> bytes_allocated_{0};
  std::atomic<int64_t> max_memory_{0};
  std::atomic<int64_t> total_bytes_allocated_{0};
  std::atomic<int64_t> num_allocations_{0};
};

/// The process-wide DecodePool every Parquet reader decodes into.
///
/// ponytail: never destroyed.  Decoded buffers can outlive main() in static
/// objects, and exit reclaims the mappings anyway.
DecodePool& decode_pool();

} // namespace MzPeak::Util

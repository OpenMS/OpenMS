/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/util/decode_pool.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstring>
#include <new>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#ifdef MADV_HUGEPAGE
#define MZPEAK_HUGE_PAGES 1
#endif
#endif

namespace MzPeak::Util {

namespace {

/// Blocks come straight from the OS: page-aligned, so 64-byte aligned, and
/// nobody's arena.  With mmap a page takes memory once it is touched (on
/// Linux, once it is written), so the part of a fresh block above the request
/// costs address space, not memory (up to huge-page rounding: a huge page is
/// allocated whole); a recycled block keeps the pages its earlier users wrote,
/// as nothing hands them back before unmap_block().  Linux does count the
/// whole mapping as committed, which matters only under strict overcommit
/// (vm.overcommit_memory=2).  VirtualAlloc(MEM_COMMIT) charges the whole
/// block, slack included, against Windows' commit limit when it is mapped;
/// there too a page takes RAM only once touched.
uint8_t* map_block(int64_t bytes)
{
  const auto len = static_cast<size_t>(bytes);
#if defined(_WIN32)
  return static_cast<uint8_t*>(
      VirtualAlloc(nullptr, len, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE));
#elif defined(MZPEAK_HUGE_PAGES)
  // Huge pages need 2 MB-aligned memory: map 2 MB more than needed and unmap
  // the unaligned ends.  With transparent huge pages set to "madvise" (the
  // usual default) only advised memory gets them.  ponytail: khugepaged may
  // later fill a partly written huge page of a long-lived block.
  const auto huge = static_cast<size_t>(DecodePool::kHugePage);
  void* raw = mmap(nullptr, len + huge, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (raw == MAP_FAILED) return nullptr;
  const auto base = reinterpret_cast<uintptr_t>(raw);
  const uintptr_t aligned = (base + huge - 1) & ~(huge - 1);
  if (aligned > base) munmap(raw, aligned - base);
  munmap(reinterpret_cast<void*>(aligned + len), base + huge - aligned);
  madvise(reinterpret_cast<void*>(aligned), len, MADV_HUGEPAGE);
  return reinterpret_cast<uint8_t*>(aligned);
#else
  void* p = mmap(nullptr, len, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  return p == MAP_FAILED ? nullptr : static_cast<uint8_t*>(p);
#endif
}

void unmap_block(void* p, int64_t bytes)
{
#ifdef _WIN32
  (void)bytes;
  VirtualFree(p, 0, MEM_RELEASE);
#else
  munmap(p, static_cast<size_t>(bytes));
#endif
}

int class_index(int64_t block)
{
  // block_size() values are 2^(e-1) + j * 2^(e-3), j in 1..4.
  const int e = std::bit_width(static_cast<uint64_t>(block - 1));
  const auto j = static_cast<int>((block - (int64_t{1} << (e - 1))) >> (e - 3));
  return 4 * (e - 1) + j - 1;
}

} // namespace

/// Written into a cached block's first bytes: the pool keeps no memory of its
/// own for the blocks it caches, so it never calls an allocator under its lock.
struct DecodePool::Node {
  Node* older;
  Node* newer; ///< across classes, oldest first
  Node* prev;
  Node* next; ///< within the class, most recent first
  int64_t size;
};

/******************************************************************************/
DecodePool::DecodePool(int64_t retain_floor)
    : system_(arrow::system_memory_pool())
    , retain_floor_(retain_floor)
{
}

/******************************************************************************/
DecodePool::~DecodePool()
{
  ReleaseUnused();
  if (region_end_ != region_next_) unmap_block(region_next_, region_end_ - region_next_);
}

/******************************************************************************/
int64_t DecodePool::block_size(int64_t size)
{
  // Four classes per power of two, never finer than 64 KB: the OS's mapping
  // granularity on Windows, and a multiple of every page size.
  const int64_t s = std::max(size, kMinBlock);
  const int e = std::bit_width(static_cast<uint64_t>(s - 1));
  const int64_t step = std::max(kMinBlock, int64_t{1} << (e - 3));
  return (s + step - 1) / step * step;
}

/******************************************************************************/
void DecodePool::roll_epoch_()
{
  // ponytail: one-second epochs of steady_clock; the peak is "recent" for one
  // to two seconds.
  const int64_t now = std::chrono::duration_cast<std::chrono::seconds>(
                          std::chrono::steady_clock::now().time_since_epoch())
                          .count();
  if (now == epoch_) return;
  last_peak_ = now == epoch_ + 1 ? peak_ : live_;
  peak_ = live_;
  epoch_ = now;
}

/******************************************************************************/
void DecodePool::unlink_(Node* n)
{
  (n->older ? n->older->newer : oldest_) = n->newer;
  (n->newer ? n->newer->older : newest_) = n->older;
  (n->prev ? n->prev->next : lists_[class_index(n->size)]) = n->next;
  if (n->next) n->next->prev = n->prev;
}

/******************************************************************************/
DecodePool::Node* DecodePool::trim_()
{
  // Checked when the live set grows as well as when a block is cached: a miss
  // that maps a fresh block while blocks of other classes sit in the cache
  // would otherwise take the footprint past the bound until the next free.
  //
  // ponytail: a quarter of the peak.  Measured at 128 threads, eviction
  // bursts free 300-500 MB against a 1.8 GB peak and blocks change class
  // under them; an eighth nearly doubled the mappings, a half bought nothing.
  Node* doomed = nullptr;
  const int64_t peak = std::max(peak_, last_peak_);
  const int64_t allowed = std::max({retain_floor_, peak / 4, peak - live_});
  while (cached_ > allowed && oldest_) {
    Node* old = oldest_;
    unlink_(old);
    cached_ -= old->size;
    old->next = doomed;
    doomed = old;
  }
  return doomed;
}

/******************************************************************************/
void DecodePool::unmap_(Node* doomed)
{
  // Outside the lock: munmap is a system call and shoots down TLBs.
  while (doomed) {
    Node* n = doomed;
    doomed = n->next;
    unmap_block(n, n->size);
  }
}

/******************************************************************************/
uint8_t* DecodePool::carve_(int64_t block)
{
  // ponytail: a new region is mapped under the lock, once per 64 MB of growth.
  if (region_end_ - region_next_ < block) {
    // What was never carved is returned: its huge page may be half-used.
    if (region_end_ != region_next_) unmap_block(region_next_, region_end_ - region_next_);
    region_next_ = map_block(kRegion);
    region_end_ = region_next_ ? region_next_ + kRegion : nullptr;
    if (!region_next_) return nullptr;
  }
  uint8_t* p = region_next_;
  region_next_ += block;
  return p;
}

/******************************************************************************/
uint8_t* DecodePool::take_(int64_t block)
{
  uint8_t* p = nullptr;
  Node* doomed = nullptr;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    roll_epoch_();
    live_ += block;
    peak_ = std::max(peak_, live_);
    if (Node* hit = lists_[class_index(block)]) {
      unlink_(hit);
      cached_ -= block;
      p = reinterpret_cast<uint8_t*>(hit);
    }
#ifdef MZPEAK_HUGE_PAGES
    else if (block < kHugePage) {
      p = carve_(block);
    }
#endif
    doomed = trim_();
  }
  unmap_(doomed);

  if (!p) p = map_block(block);
  if (!p) {
    // Out of address space or memory: what the cache holds may be the room.
    ReleaseUnused();
    p = map_block(block);
  }
  if (!p) {
    std::lock_guard<std::mutex> lock(mutex_);
    live_ -= block;
  }
  return p;
}

/******************************************************************************/
void DecodePool::give_(uint8_t* p, int64_t block)
{
  Node* doomed = nullptr;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    roll_epoch_();
    live_ -= block;

    Node*& head = lists_[class_index(block)];
    Node* n = new (p) Node{newest_, nullptr, nullptr, head, block};
    (newest_ ? newest_->newer : oldest_) = n;
    newest_ = n;
    if (head) head->prev = n;
    head = n;
    cached_ += block;
    doomed = trim_();
  }
  unmap_(doomed);
}

/******************************************************************************/
void DecodePool::ReleaseUnused()
{
  Node* doomed = nullptr;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    doomed = oldest_;
    oldest_ = newest_ = nullptr;
    std::fill(std::begin(lists_), std::end(lists_), nullptr);
    cached_ = 0;
  }
  while (doomed) {
    Node* n = doomed;
    doomed = n->newer;
    unmap_block(n, n->size);
  }
}

/******************************************************************************/
arrow::Status DecodePool::Allocate(int64_t size, int64_t alignment, uint8_t** out)
{
  if (size < 0) return arrow::Status::Invalid("negative malloc size");
  if (!block_backed_(size, alignment)) {
    ARROW_RETURN_NOT_OK(system_->Allocate(size, alignment, out));
  } else {
    uint8_t* p = take_(block_size(size));
    if (!p) return arrow::Status::OutOfMemory("DecodePool: cannot map ", size, " bytes");
    *out = p;
  }
  did_allocate_(size);
  return arrow::Status::OK();
}

/******************************************************************************/
arrow::Status DecodePool::Reallocate(int64_t old_size,
                                     int64_t new_size,
                                     int64_t alignment,
                                     uint8_t** ptr)
{
  if (new_size < 0) return arrow::Status::Invalid("negative realloc size");
  const bool was_block = block_backed_(old_size, alignment);
  const bool is_block = block_backed_(new_size, alignment);

  if (!was_block && !is_block) {
    ARROW_RETURN_NOT_OK(system_->Reallocate(old_size, new_size, alignment, ptr));
  } else if (was_block && is_block && block_size(old_size) == block_size(new_size)) {
    // Same class: the block already has the room.  A grow-then-shrink-to-fit
    // column buffer (Parquet reserves double, then trims) mostly stays here.
  } else {
    // Another class: move, so a block's size always follows from the size
    // Arrow passes back.  That keeps the pool free of any per-block lookup,
    // and a big block never lingers behind a small array.
    uint8_t* fresh = nullptr;
    if (is_block) {
      fresh = take_(block_size(new_size));
      if (!fresh) return arrow::Status::OutOfMemory("DecodePool: cannot map ", new_size, " bytes");
    } else {
      ARROW_RETURN_NOT_OK(system_->Allocate(new_size, alignment, &fresh));
    }
    std::memcpy(fresh, *ptr, static_cast<std::size_t>(std::min(old_size, new_size)));
    if (was_block) {
      give_(*ptr, block_size(old_size));
    } else {
      system_->Free(*ptr, old_size, alignment);
    }
    *ptr = fresh;
  }

  if (new_size > old_size) {
    did_allocate_(new_size - old_size);
  } else {
    did_free_(old_size - new_size);
  }
  return arrow::Status::OK();
}

/******************************************************************************/
void DecodePool::Free(uint8_t* buffer, int64_t size, int64_t alignment)
{
  if (block_backed_(size, alignment)) {
    give_(buffer, block_size(size));
  } else {
    system_->Free(buffer, size, alignment);
  }
  did_free_(size);
}

/******************************************************************************/
void DecodePool::did_allocate_(int64_t bytes)
{
  // As Arrow's own pools count: a reallocation that grows is one allocation
  // of the difference.
  const int64_t now = bytes_allocated_.fetch_add(bytes, std::memory_order_relaxed) + bytes;
  total_bytes_allocated_.fetch_add(bytes, std::memory_order_relaxed);
  num_allocations_.fetch_add(1, std::memory_order_relaxed);
  int64_t peak = max_memory_.load(std::memory_order_relaxed);
  while (peak < now &&
         !max_memory_.compare_exchange_weak(peak, now, std::memory_order_relaxed)) {
  }
}

/******************************************************************************/
void DecodePool::did_free_(int64_t bytes)
{
  bytes_allocated_.fetch_sub(bytes, std::memory_order_relaxed);
}

/******************************************************************************/
int64_t DecodePool::bytes_allocated() const { return bytes_allocated_.load(); }
int64_t DecodePool::max_memory() const { return max_memory_.load(); }
int64_t DecodePool::total_bytes_allocated() const { return total_bytes_allocated_.load(); }
int64_t DecodePool::num_allocations() const { return num_allocations_.load(); }
std::string DecodePool::backend_name() const { return "mzpeak-decode(os-blocks+system)"; }

/******************************************************************************/
int64_t DecodePool::bytes_mapped() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return live_ + cached_;
}

/******************************************************************************/
int64_t DecodePool::bytes_cached() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return cached_;
}

/******************************************************************************/
DecodePool& decode_pool()
{
  static DecodePool* pool = new DecodePool;
  return *pool;
}

} // namespace MzPeak::Util

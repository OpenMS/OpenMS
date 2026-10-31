/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <tuple>
#include <utility>
#include <vector>

#ifdef __linux__
#include <sys/mman.h>
#include <unistd.h>
#endif

namespace MzPeak::Util {

/**
 * A sorted, contiguous entity-index -> value map.
 *
 * WHY NOT std::map.  The per-spectrum metadata is one entry per spectrum and
 * the entries are read once and then only looked up, so a node-based tree pays
 * for an insertion order it never needs: 32 bytes of pointers and colour per
 * entry, and -- the part that actually costs -- one malloc per entry.  A run
 * with 300k spectra means 300k separate small allocations whose rounding and
 * fragmentation the process never gets back (macOS does not return freed small
 * blocks to the OS).  One vector is one allocation.
 *
 * Lookup is better too, not merely no worse.  Entity indices are dense and
 * 0-based in every archive this reader has seen, so find() hits `data_[key]`
 * directly -- O(1), one cache line -- and falls back to a binary search when a
 * file does not oblige.  The fallback is what makes the fast path safe to try.
 *
 * USE.  Fill with append() in file order, call sort() ONCE, then look up.
 * There is deliberately no map-like operator[] insert: an insert into a sorted
 * vector is O(n), and offering one would hide that behind familiar syntax.
 * Mutating a value found by find() is fine and is how the facet passes attach
 * precursors and scans.
 */
template <typename T> class IndexMap final {
public:
  using value_type = std::pair<uint64_t, T>;
  using container_type = std::vector<value_type>;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

  iterator begin() { return data_.begin(); }
  iterator end() { return data_.end(); }
  const_iterator begin() const { return data_.begin(); }
  const_iterator end() const { return data_.end(); }
  /// Reserve room for @p n entries -- on Linux, backed by transparent huge
  /// pages when that is large.
  ///
  /// WHY: the map is filled once, front to back, right after this, and on
  /// 4 KB pages that fill is mostly page faults: 76k of them for a
  /// 717,924-spectrum run (310 MB), which made PASS 1 of the metadata read
  /// 0.30 s instead of 0.16 s on kim.  2 MB pages take 512x fewer, and
  /// unmapping the map at the end gets the same factor.
  ///
  /// ponytail: advice, never a requirement -- a kernel that declines it (THP
  /// "never", or no huge page free) keeps 4 KB pages and loses nothing, and
  /// below 32 MB the faults are not worth a syscall.
  void reserve(std::size_t n)
  {
    data_.reserve(n);
#ifdef __linux__
    const std::size_t bytes = data_.capacity() * sizeof(value_type);
    if (bytes < (std::size_t(32) << 20)) return;
    const auto page = static_cast<std::uintptr_t>(sysconf(_SC_PAGESIZE));
    const auto begin = reinterpret_cast<std::uintptr_t>(data_.data());
    const std::uintptr_t first = (begin + page - 1) & ~(page - 1);
    const std::uintptr_t last = (begin + bytes) & ~(page - 1);
    if (last > first) madvise(reinterpret_cast<void*>(first), last - first, MADV_HUGEPAGE);
#endif
  }

  /// Append in file order; sort() must follow before any find().
  void append(uint64_t key, T value) { data_.emplace_back(key, std::move(value)); }

  /// Append a default value in file order and return it, to be filled in
  /// place: building a large value on the side and moving it in costs a
  /// copy of it per entry.
  T& append(uint64_t key)
  {
    return data_.emplace_back(std::piecewise_construct, std::forward_as_tuple(key),
                              std::forward_as_tuple()).second;
  }

  /// Order by key, keeping the FIRST entry for a duplicated key.
  ///
  /// Checked before sorted: entities are written in ascending index order, so
  /// the sort is normally not needed at all -- and skipping it matters for
  /// memory, not just time, because stable_sort allocates a temporary the size
  /// of the whole map.  On a 7,534-spectrum run that temporary was 3.3 MB of
  /// peak that the steady state never shows.
  ///
  /// stable_sort when it IS needed, not sort: a file that repeats an index
  /// would otherwise resolve to whichever duplicate the sort happened to leave
  /// first, which is a different spectrum from run to run.  Duplicates are not
  /// expected, but "not expected" is not a reason to be nondeterministic.
  void sort()
  {
    const auto by_key = [](const value_type& a, const value_type& b) {
      return a.first < b.first;
    };
    if (!std::is_sorted(data_.begin(), data_.end(), by_key))
      std::stable_sort(data_.begin(), data_.end(), by_key);

    data_.erase(std::unique(data_.begin(), data_.end(),
                            [](const value_type& a, const value_type& b) {
                              return a.first == b.first;
                            }),
                data_.end());
  }

  iterator find(uint64_t key)
  {
    const const_iterator it = static_cast<const IndexMap&>(*this).find(key);
    return data_.begin() + (it - data_.cbegin());
  }

  const_iterator find(uint64_t key) const
  {
    // Dense 0-based indices: the key IS the slot.  Checked, not assumed.
    if (key < data_.size() && data_[static_cast<std::size_t>(key)].first == key)
      return data_.cbegin() + static_cast<std::ptrdiff_t>(key);

    const const_iterator it = std::lower_bound(
        data_.cbegin(), data_.cend(), key,
        [](const value_type& e, uint64_t k) { return e.first < k; });
    return (it != data_.cend() && it->first == key) ? it : data_.cend();
  }

private:
  container_type data_;
};

} // namespace MzPeak::Util

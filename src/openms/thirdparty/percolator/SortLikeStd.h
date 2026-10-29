/*
 * OpenMS addition (not part of upstream Percolator): sorts like std::sort, with the tasks of
 * the enclosing OpenMP parallel region.
 *
 * libstdc++'s std::sort (bits/stl_algo.h, std::__sort) is an introsort: __introsort_loop
 * partitions a range around the median of three elements, recurses into the right part, loops
 * on the left part and stops at ranges of at most _S_threshold = 16 elements (heap sort once
 * 2 * floor(log2(n)) partitionings are nested); one insertion sort over the whole range
 * (__final_insertion_sort) finishes the job. sortLikeStdSort makes the same comparisons and
 * swaps, so it leaves the same permutation - also among elements that compare equal - but
 * hands the right parts to OpenMP tasks:
 *  - a partitioning step reads and writes only its own range (the scans of
 *    __unguarded_partition are stopped by elements of the range), and the two parts are
 *    disjoint, so the order in which they are processed does not matter;
 *  - after the partitioning the ranges of <= 16 elements ("leaves") are in ascending order
 *    relative to each other, so the final insertion sort never moves an element across a leaf
 *    boundary: it is an insertion sort of each leaf, which can be done as soon as the leaf is
 *    known. (__insertion_sort, used for the first 16 elements, and __unguarded_linear_insert
 *    differ in their comparisons, not in where an element ends up.)
 * Both points need what std::sort needs anyway: a comparison that is a strict weak ordering
 * and has no state. The partitionings of large ranges (the first ones, where the scans of one
 * partitioning are the critical path) make the swaps of __unguarded_partition with all tasks
 * (partitionInParallel); partitionings of disjoint ranges may run so at the same time.
 * The same algorithm sorts the peptides of OpenMS' FragmentIndex (FragmentIndex::sortPeptides_,
 * whose test compares it with std::sort).
 *
 * Only with libstdc++ and OpenMP, and only inside an active parallel region (the tasks run on
 * its threads); otherwise, and for small ranges, it calls std::sort.
 */
#ifndef PERCOLATOR_SORT_LIKE_STD_H
#define PERCOLATOR_SORT_LIKE_STD_H

#include <algorithm>
#include <bit>
#include <cstddef>
#include <functional>
#include <iterator>
#include <new>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

namespace SortLikeStdDetail {

constexpr std::ptrdiff_t kLeafSize = 16;      // libstdc++'s _S_threshold
constexpr std::ptrdiff_t kMinTaskSize = 2048;  // right parts at least this large become tasks
// ranges of more elements are partitioned by partitionInParallel: only large sets (at smaller
// sizes, the extra pass of partitionInParallel competed with the sorts of the other folds)
constexpr std::ptrdiff_t kMinParallelPartition = std::ptrdiff_t(1) << 18;
constexpr std::ptrdiff_t kChunksPerPartition = 64;  // chunks of partitionInParallel's counting pass
constexpr std::ptrdiff_t kMinChunk = 2048;          // elements a task of partitionInParallel scans at least
constexpr std::ptrdiff_t kSwapBlocksPerPartition = 16;
constexpr std::ptrdiff_t kMinSwapBlock = 1024;

#if defined(__GLIBCXX__) && defined(_OPENMP)
// std::__unguarded_partition(first + 1, last, first) with the tasks of the enclosing parallel
// region, for large ranges, where its sequential scans are the critical path of the sort (the
// idea of FragmentIndex::sortPeptides_). The scans swap the k-th element from the left that is
// not less than the pivot (a_k) with the k-th element from the right that the pivot is not less
// than (b_k), for k = 0, 1, ... while a_k < b_k, and stop at K, the first k with a_k >= b_k: the
// scan from the left stops at a_K or, if that lies beyond, at b_{K-1}, which holds a_{K-1}'s
// element by then. (Each scan only passes elements it would not stop at, which no swap has
// touched, up to the place of the last swap.) So the same swaps are found from where both kinds
// of elements are: K is the largest k for which some position m has at least k of the first
// kind left of it and k of the second kind from m on. The kinds are compared once, in chunks of
// the range (in parallel), and noted per element, so that a_k and b_k are found without
// comparing again. Returns where the scan from the left stops, or nullptr (nothing changed) if
// its bookkeeping cannot be allocated.
template <typename T, typename Less>
T* partitionInParallel(T* first, T* last, Less less, std::ptrdiff_t chunk) {
  T* const begin = first + 1;  // the pivot *first is not moved: the scans cover [first + 1, last)
  const std::ptrdiff_t n = last - begin;
  const std::ptrdiff_t numChunks = n / chunk + (n % chunk != 0);
  const auto chunkBegin = [n, chunk, numChunks](std::ptrdiff_t c) { return c < numChunks ? c * chunk : n; };
  // kinds[i] & 1: the scan from the left stops at begin[i] (!less(begin[i], pivot));
  // kinds[i] & 2: the scan from the right stops there (!less(pivot, begin[i]))
  std::vector<unsigned char> kinds;
  // left[c]: elements of the first kind in chunks [0, c); right[c]: of the second kind in
  // chunks [c, numChunks)
  std::vector<std::ptrdiff_t> left, right;
  try {
    kinds.resize(static_cast<std::size_t>(n));
    left.assign(numChunks + 1, 0);
    right.assign(numChunks + 1, 0);
  } catch (const std::bad_alloc&) {
    return nullptr;
  }
#pragma omp taskloop grainsize(1) default(none) shared(kinds, left, right, less, chunkBegin, numChunks, first, begin)
  for (std::ptrdiff_t c = 0; c < numChunks; ++c) {
    const T& pivot = *first;
    std::ptrdiff_t l = 0, r = 0;
    for (std::ptrdiff_t i = chunkBegin(c); i != chunkBegin(c + 1); ++i) {
      const unsigned char kind = static_cast<unsigned char>(!less(begin[i], pivot)) |
                                 static_cast<unsigned char>(!less(pivot, begin[i]) << 1);
      kinds[i] = kind;
      l += kind & 1;
      r += kind >> 1;
    }
    left[c + 1] = l;
    right[c] = r;
  }
  for (std::ptrdiff_t c = 0; c < numChunks; ++c) left[c + 1] += left[c];
  for (std::ptrdiff_t c = numChunks; c-- > 0;) right[c] += right[c + 1];

  // K = max over m of min(#first kind in [begin, m), #second kind in [m, last)); the first
  // count grows with m, the second shrinks, so the maximum lies where they cross: in the chunk
  // where the first passes the second, or at one of its ends
  std::ptrdiff_t cross = 0;  // first chunk boundary with left >= right
  while (cross < numChunks && left[cross] < right[cross]) ++cross;
  std::ptrdiff_t K = std::min(left[cross], right[cross]);
  if (cross > 0) {
    std::ptrdiff_t l = left[cross - 1], r = right[cross - 1];
    K = std::max(K, std::min(l, r));
    for (std::ptrdiff_t i = chunkBegin(cross - 1); i != chunkBegin(cross); ++i) {
      l += kinds[i] & 1;
      r -= kinds[i] >> 1;
      K = std::max(K, std::min(l, r));
    }
  }

  // positions (in begin[]) of a_k and b_k (k < number of such elements; otherwise n, and -1:
  // no such element)
  const auto aAt = [&](std::ptrdiff_t k) -> std::ptrdiff_t {
    if (k >= left[numChunks]) return n;
    const std::ptrdiff_t c = std::upper_bound(left.begin(), left.end(), k) - left.begin() - 1;
    std::ptrdiff_t seen = left[c];
    for (std::ptrdiff_t i = chunkBegin(c);; ++i) {
      if ((kinds[i] & 1) && seen++ == k) return i;
    }
  };
  const auto bAt = [&](std::ptrdiff_t k) -> std::ptrdiff_t {
    if (k >= right[0]) return -1;
    // the chunk c with right[c + 1] <= k < right[c]
    const std::ptrdiff_t c = std::lower_bound(right.begin(), right.end(), k,
                                              std::greater<std::ptrdiff_t>()) - right.begin() - 1;
    std::ptrdiff_t seen = right[c + 1];
    for (std::ptrdiff_t i = chunkBegin(c + 1) - 1;; --i) {
      if ((kinds[i] >> 1) && seen++ == k) return i;
    }
  };
  T* const cut = begin + ((K == 0) ? aAt(0) : std::min(aAt(K), bAt(K - 1)));

  // the swaps, in blocks of consecutive k
  const std::ptrdiff_t block = std::max<std::ptrdiff_t>(
      kMinSwapBlock, (K + kSwapBlocksPerPartition - 1) / kSwapBlocksPerPartition);
  const std::ptrdiff_t numBlocks = (K + block - 1) / block;
#pragma omp taskloop grainsize(1) default(none) shared(kinds, aAt, bAt, numBlocks, block, K, begin)
  for (std::ptrdiff_t j = 0; j < numBlocks; ++j) {
    std::ptrdiff_t a = aAt(j * block);
    std::ptrdiff_t b = bAt(j * block);
    for (std::ptrdiff_t k = j * block; k < std::min(K, (j + 1) * block); ++k) {
      if (k > j * block) {
        while (!(kinds[++a] & 1)) {}
        while (!(kinds[--b] >> 1)) {}
      }
      std::iter_swap(begin + a, begin + b);
    }
  }
  return cut;
}

template <typename T, typename Less>
void introsortLoop(T* first, T* last, int depth_limit, Less less,
                   std::ptrdiff_t minParallelPartition) {
  while (last - first > kLeafSize) {
    if (depth_limit == 0) {
      // std::__partial_sort(first, last, last, comp) == __make_heap + __sort_heap
      std::make_heap(first, last, less);
      std::sort_heap(first, last, less);
      return;
    }
    --depth_limit;
    // std::__unguarded_partition_pivot: __move_median_to_first(first, first + 1, mid, last - 1) ...
    T* const a = first + 1;
    T* const b = first + (last - first) / 2;
    T* const c = last - 1;
    if (less(*a, *b)) {
      if (less(*b, *c)) std::iter_swap(first, b);
      else if (less(*a, *c)) std::iter_swap(first, c);
      else std::iter_swap(first, a);
    } else if (less(*a, *c)) std::iter_swap(first, a);
    else if (less(*b, *c)) std::iter_swap(first, c);
    else std::iter_swap(first, b);
    // ... and __unguarded_partition(first + 1, last, pivot = first)
    T* cut = nullptr;
    if (last - first > minParallelPartition) {
      // chunks: at least kMinChunk elements, at most about kChunksPerPartition
      cut = partitionInParallel(first, last, less,
                                std::max(kMinChunk, (last - first) / kChunksPerPartition + 1));
    }
    if (cut == nullptr) {
      cut = first + 1;
      T* hi = last;
      while (true) {
        while (less(*cut, *first)) ++cut;
        --hi;
        while (less(*first, *hi)) --hi;
        if (!(cut < hi)) break;
        std::iter_swap(cut, hi);
        ++cut;
      }
    }
    // right part [cut, last): recursion in std::sort, a task here if it is large enough
    if (last - cut >= kMinTaskSize) {
#pragma omp task default(none) firstprivate(cut, last, depth_limit, less, minParallelPartition)
      introsortLoop(cut, last, depth_limit, less, minParallelPartition);
    } else {
      introsortLoop(cut, last, depth_limit, less, minParallelPartition);
    }
    last = cut;  // left part [first, cut): next iteration
  }
  // leaf: the part of __final_insertion_sort that concerns [first, last)
  for (T* i = first + 1; i < last; ++i) {
    T value = std::move(*i);
    T* pos = i;
    for (; pos != first && less(value, *(pos - 1)); --pos) *pos = std::move(*(pos - 1));
    *pos = std::move(value);
  }
}
#endif

// sortLikeStdSort, with ranges of more than minParallelPartition elements partitioned by
// partitionInParallel (for tests)
template <typename T, typename Less>
void sortLikeStdSort(T* first, T* last, Less less, std::ptrdiff_t minParallelPartition) {
#if defined(__GLIBCXX__) && defined(_OPENMP)
  const std::ptrdiff_t n = last - first;
  if (n > 2 * kMinTaskSize && omp_in_parallel()) {
    const int depth_limit = 2 * (static_cast<int>(std::bit_width(static_cast<std::size_t>(n))) - 1);  // std::__lg(n) * 2
#pragma omp taskgroup
    introsortLoop(first, last, depth_limit, less, minParallelPartition);
    return;
  }
#else
  (void)minParallelPartition;
#endif
  std::sort(first, last, less);
}

}  // namespace SortLikeStdDetail

// Sorts [first, last) of a contiguous range exactly like std::sort(first, last, less).
// The comparison must be a strict weak ordering without state (as for std::sort) and must
// not throw; T's moves and swaps must not throw.
template <typename T, typename Less>
void sortLikeStdSort(T* first, T* last, Less less) {
  SortLikeStdDetail::sortLikeStdSort(first, last, less, SortLikeStdDetail::kMinParallelPartition);
}

// The whole vector (or other contiguous container).
template <typename Container, typename Less>
void sortLikeStdSort(Container& v, Less less) {
  if (v.empty()) return;
  sortLikeStdSort(v.data(), v.data() + v.size(), less);
}

}}}  // namespace OpenMS::Internal::Percolator

#endif  // PERCOLATOR_SORT_LIKE_STD_H

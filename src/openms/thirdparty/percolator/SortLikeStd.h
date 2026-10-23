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
 * and has no state. The same algorithm sorts the peptides of OpenMS' FragmentIndex
 * (FragmentIndex::sortPeptides_, whose test compares it with std::sort).
 *
 * Only with libstdc++ and OpenMP, and only inside an active parallel region (the tasks run on
 * its threads); otherwise, and for small ranges, it calls std::sort.
 */
#ifndef PERCOLATOR_SORT_LIKE_STD_H
#define PERCOLATOR_SORT_LIKE_STD_H

#include <algorithm>
#include <bit>
#include <cstddef>
#include <iterator>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

namespace SortLikeStdDetail {

constexpr std::ptrdiff_t kLeafSize = 16;      // libstdc++'s _S_threshold
constexpr std::ptrdiff_t kMinTaskSize = 2048;  // right parts at least this large become tasks

#if defined(__GLIBCXX__) && defined(_OPENMP)
template <typename T, typename Less>
void introsortLoop(T* first, T* last, int depth_limit, Less less) {
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
    T* cut = first + 1;
    T* hi = last;
    while (true) {
      while (less(*cut, *first)) ++cut;
      --hi;
      while (less(*first, *hi)) --hi;
      if (!(cut < hi)) break;
      std::iter_swap(cut, hi);
      ++cut;
    }
    // right part [cut, last): recursion in std::sort, a task here if it is large enough
    if (last - cut >= kMinTaskSize) {
#pragma omp task default(none) firstprivate(cut, last, depth_limit, less)
      introsortLoop(cut, last, depth_limit, less);
    } else {
      introsortLoop(cut, last, depth_limit, less);
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

}  // namespace SortLikeStdDetail

// Sorts [first, last) of a contiguous range exactly like std::sort(first, last, less).
// The comparison must be a strict weak ordering without state (as for std::sort) and must
// not throw; T's moves and swaps must not throw.
template <typename T, typename Less>
void sortLikeStdSort(T* first, T* last, Less less) {
#if defined(__GLIBCXX__) && defined(_OPENMP)
  const std::ptrdiff_t n = last - first;
  if (n > 2 * SortLikeStdDetail::kMinTaskSize && omp_in_parallel()) {
    const int depth_limit = 2 * (static_cast<int>(std::bit_width(static_cast<std::size_t>(n))) - 1);  // std::__lg(n) * 2
#pragma omp taskgroup
    SortLikeStdDetail::introsortLoop(first, last, depth_limit, less);
    return;
  }
#endif
  std::sort(first, last, less);
}

// The whole vector (or other contiguous container).
template <typename Container, typename Less>
void sortLikeStdSort(Container& v, Less less) {
  if (v.empty()) return;
  sortLikeStdSort(v.data(), v.data() + v.size(), less);
}

}}}  // namespace OpenMS::Internal::Percolator

#endif  // PERCOLATOR_SORT_LIKE_STD_H

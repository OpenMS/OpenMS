/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
#ifndef PARALLEL_SORT_H_
#define PARALLEL_SORT_H_

#include <algorithm>
#include <cstddef>
#include <exception>
#include <iterator>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "OpenMPThreads.h"  // OpenMS: maxThreads(), threadsFor()

namespace OpenMS { namespace Internal { namespace Percolator {


/*
 * parallelSort(first, last, comp) sorts like std::sort, on the OpenMP
 * threads. For a comparator that is a total order, as the comparators of
 * all PSM-sized sorts are, the result is the same as std::sort's on any
 * number of threads. Elements that compare equal must then be
 * indistinguishable: -0.0 and 0.0, for example, compare equal but print
 * differently, and their order would depend on the number of threads.
 *
 * The range is split into one part per thread, like a quicksort: each split
 * partitions a range around the element at the wanted quantile of an evenly
 * spaced sample, with every thread partitioning a chunk and then swapping
 * the misplaced elements across the split point. The parts are then sorted
 * with std::sort. All steps run in parallel and work in place, so the sort
 * needs no memory besides the range (GCC's parallel mode copied it).
 *
 * Only "parallel for" over int indices is used, so that it also works with
 * OpenMP 2.0 (MSVC). Without OpenMP, inside a parallel region, with one
 * thread or for short ranges it is std::sort. An exception thrown by the
 * comparator on any thread is passed on to the caller (as from std::sort,
 * the range is then left in an unspecified state).
 */
#ifdef _OPENMP
namespace parallel_sort_detail {

// below this many elements a parallel region costs more than it saves
const std::ptrdiff_t kMinParallelSize = 1 << 15;
// largest sample for choosing a split element; with s elements the parts
// deviate from their wanted size by about 1 / (2 sqrt(s)) of the split range
const std::ptrdiff_t kMaxSampleSize = 1 << 14;

// x * num / den without overflow for ranges of up to 2^31 elements also where
// ptrdiff_t has 32 bits
inline std::ptrdiff_t scaled(std::ptrdiff_t x, long long num, long long den) {
  return static_cast<std::ptrdiff_t>(static_cast<long long>(x) * num / den);
}

// Runs f(i) for i in [0, n) on the OpenMP threads and rethrows the first
// exception after the loop.
template <class F>
void parallelFor(int n, const F& f) {
  std::exception_ptr error;
#pragma omp parallel for schedule(static) num_threads(threadsFor(n))
  for (int i = 0; i < n; ++i) {
    try {
      f(i);
    } catch (...) {
#pragma omp critical(parallel_sort_error)
      if (!error) error = std::current_exception();
    }
  }
  if (error) std::rethrow_exception(error);
}

// Consecutive pieces [begin, end) of the range, numbered on through them
// by the running total "before" of the preceding pieces' lengths.
struct Segment {
  std::ptrdiff_t begin, end, before;
};

// Swaps the k-th element of the pieces a with the k-th of the pieces b, for
// k in [from, to). Both lists cover at least to elements.
template <class RandomIt>
void swapSegments(RandomIt first, const std::vector<Segment>& a,
                  const std::vector<Segment>& b, std::ptrdiff_t from,
                  std::ptrdiff_t to) {
  auto findSegment = [from](const std::vector<Segment>& s) {
    std::size_t i = 0;
    while (s[i].before + (s[i].end - s[i].begin) <= from) ++i;
    return i;
  };
  std::size_t ia = findSegment(a), ib = findSegment(b);
  std::ptrdiff_t pa = a[ia].begin + (from - a[ia].before);
  std::ptrdiff_t pb = b[ib].begin + (from - b[ib].before);
  for (std::ptrdiff_t k = from; k < to; ++k) {
    while (pa == a[ia].end) pa = a[++ia].begin;
    while (pb == b[ib].end) pb = b[++ib].begin;
    std::iter_swap(first + pa, first + pb);
    ++pa;
    ++pb;
  }
}

// Partitions [first + lo, first + hi) on p threads into the elements that
// compare less than a split element and the others, and returns the offset
// of the first of the others. About (hi - lo) * num / den elements are less
// than the split element.
template <class RandomIt, class Compare>
std::ptrdiff_t parallelSplit(RandomIt first, std::ptrdiff_t lo,
                             std::ptrdiff_t hi, int num, int den, int p,
                             Compare comp) {
  typedef typename std::iterator_traits<RandomIt>::value_type Value;
  const std::ptrdiff_t len = hi - lo;
  if (len < kMinParallelSize) {  // a short range is split exactly
    const std::ptrdiff_t split = lo + scaled(len, num, den);
    std::nth_element(first + lo, first + split, first + hi, comp);
    return split;
  }
  // the split element: the wanted quantile of an evenly spaced sample
  const std::ptrdiff_t numSample = std::min(len / 32, kMaxSampleSize);
  std::vector<Value> sample;
  sample.reserve(static_cast<std::size_t>(numSample));
  for (std::ptrdiff_t k = 0; k < numSample; ++k) {
    sample.push_back(first[lo + scaled(len, 2 * k + 1, 2 * numSample)]);
  }
  const std::ptrdiff_t rank = scaled(numSample, num, den);
  std::nth_element(sample.begin(), sample.begin() + rank, sample.end(), comp);
  const Value pivot = sample[static_cast<std::size_t>(rank)];
  std::vector<Value>().swap(sample);

  // partition chunk c, [lo + len * c / p, lo + len * (c + 1) / p), with the
  // smaller elements first
  std::vector<std::ptrdiff_t> chunk(static_cast<std::size_t>(p) + 1),
      numSmall(static_cast<std::size_t>(p));
  for (int c = 0; c <= p; ++c) {
    chunk[static_cast<std::size_t>(c)] = lo + scaled(len, c, p);
  }
  parallelFor(p, [&](int c) {
    const std::size_t i = static_cast<std::size_t>(c);
    RandomIt begin = first + chunk[i];
    RandomIt mid = std::partition(
        begin, first + chunk[i + 1],
        [&](const Value& x) { return comp(x, pivot); });
    numSmall[i] = mid - begin;
  });
  std::ptrdiff_t split = lo;
  for (int c = 0; c < p; ++c) split += numSmall[static_cast<std::size_t>(c)];

  // The larger elements before split and the smaller ones after it are
  // equally many; swap them pairwise.
  std::vector<Segment> large, small;
  std::ptrdiff_t numLarge = 0, numMisplacedSmall = 0;
  for (int c = 0; c < p; ++c) {
    const std::size_t i = static_cast<std::size_t>(c);
    const std::ptrdiff_t mid = chunk[i] + numSmall[i];
    const std::ptrdiff_t largeEnd = std::min(chunk[i + 1], split);
    if (mid < largeEnd) {
      Segment s = {mid, largeEnd, numLarge};
      large.push_back(s);
      numLarge += largeEnd - mid;
    }
    const std::ptrdiff_t smallBegin = std::max(chunk[i], split);
    if (smallBegin < mid) {
      Segment s = {smallBegin, mid, numMisplacedSmall};
      small.push_back(s);
      numMisplacedSmall += mid - smallBegin;
    }
  }
  if (numLarge > 0) {
    parallelFor(p, [&](int t) {
      const std::ptrdiff_t from = scaled(numLarge, t, p);
      const std::ptrdiff_t to = scaled(numLarge, t + 1, p);
      if (from < to) swapSegments(first, large, small, from, to);
    });
  }
  return split;
}

}  // namespace parallel_sort_detail
#endif

template <class RandomIt, class Compare>
void parallelSort(RandomIt first, RandomIt last, Compare comp) {
#ifdef _OPENMP
  using parallel_sort_detail::parallelFor;
  const std::ptrdiff_t n = last - first;
  const int p = inParallel() ? 1 : maxThreads();  // OpenMS: the cap of the calling thread (OpenMPThreads.h)
  if (p < 2 || n < parallel_sort_detail::kMinParallelSize) {
    std::sort(first, last, comp);
    return;
  }
  // Part i is [first + bound[i], first + bound[i + 1]). Groups of parts
  // [lo, hi) are split in halves, one level after the other, each split by
  // all threads.
  std::vector<std::ptrdiff_t> bound(static_cast<std::size_t>(p) + 1, 0);
  bound[static_cast<std::size_t>(p)] = n;
  std::vector<std::pair<int, int> > groups(1, std::make_pair(0, p)), next;
  while (!groups.empty()) {
    next.clear();
    for (const std::pair<int, int>& group : groups) {
      const int lo = group.first, hi = group.second, mid = (lo + hi) / 2;
      bound[static_cast<std::size_t>(mid)] =
          parallel_sort_detail::parallelSplit(
              first, bound[static_cast<std::size_t>(lo)],
              bound[static_cast<std::size_t>(hi)], mid - lo, hi - lo, p, comp);
      if (mid - lo > 1) next.push_back(std::make_pair(lo, mid));
      if (hi - mid > 1) next.push_back(std::make_pair(mid, hi));
    }
    groups.swap(next);
  }
  parallelFor(p, [&](int i) {
    std::sort(first + bound[static_cast<std::size_t>(i)],
              first + bound[static_cast<std::size_t>(i) + 1], comp);
  });
#else
  std::sort(first, last, comp);
#endif
}

}}}  // namespace OpenMS::Internal::Percolator

#endif /* PARALLEL_SORT_H_ */

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// NOT part of the upstream Percolator source tree. See header for provenance
// and re-sync instructions.

#include "InProcessSampler.h"
#include "PseudoRandom.h"

#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <utility>

namespace OpenMS { namespace Internal { namespace Percolator {
using namespace std;

std::vector<std::size_t> sampleTrainingRowIndices(
    const std::vector<int>& scan_numbers,
    const std::vector<double>& exp_masses,
    std::size_t max_train)
{
  const std::size_t n_rows = scan_numbers.size();
  std::vector<std::size_t> out;

  if (max_train == 0 || max_train >= n_rows)
  {
    out.resize(n_rows);
    std::iota(out.begin(), out.end(), std::size_t{0});
    return out;
  }

  // ScanId = (scan, exp_mass). Same-ScanId rows share a single random
  // priority (upstream SetHandler.cpp:249-258).
  using ScanId = std::pair<int, double>;
  std::map<ScanId, std::size_t> scan_pri;

  // Max-heap by priority. Default std::priority_queue is max-first on the
  // pair's first element, which matches the upstream semantics.
  std::priority_queue<std::pair<std::size_t, std::size_t>> heap;
  std::size_t upper_limit = std::numeric_limits<std::size_t>::max();

  for (std::size_t i = 0; i < n_rows; ++i)
  {
    const ScanId sid{scan_numbers[i], exp_masses[i]};
    std::size_t pri;
    auto it = scan_pri.find(sid);
    if (it != scan_pri.end())
    {
      pri = it->second;
    }
    else
    {
      pri = PseudoRandom::lcg_rand();
      scan_pri.emplace(sid, pri);
    }

    // Upstream admission test: push iff room or below upper_limit.
    if (heap.size() < max_train || pri < upper_limit)
    {
      heap.emplace(pri, i);
      if (heap.size() > max_train)
      {
        // Upstream records the DROPPED priority as upper_limit (not the new
        // heap top). See SetHandler.cpp:267-272.
        upper_limit = heap.top().first;
        heap.pop();
      }
    }
  }

  out.reserve(heap.size());
  while (!heap.empty())
  {
    out.push_back(heap.top().second);
    heap.pop();
  }
  std::sort(out.begin(), out.end());
  return out;
}

}}}  // namespace OpenMS::Internal::Percolator

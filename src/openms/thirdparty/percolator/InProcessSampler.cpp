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
#include <stdexcept>
#include <unordered_map>
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

namespace
{
  // SplitMix64's finalizer: a bijective mix of 64 bits
  std::uint64_t mix64(std::uint64_t x)
  {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
  }
}

OutOfFoldSample drawOutOfFoldSample(const std::vector<int>& spec_files,
                                   const std::vector<int>& scans,
                                   std::size_t max_train,
                                   std::uint64_t seed,
                                   unsigned int num_folds)
{
  if (spec_files.size() != scans.size() || max_train == 0 || num_folds == 0 || num_folds > 255)
  {
    throw std::invalid_argument("drawOutOfFoldSample: inconsistent arguments");
  }
  const std::size_t n_rows = scans.size();
  struct Spectrum
  {
    std::uint64_t key;       // spec_file (high 32 bits), scan (low 32 bits)
    std::uint64_t priority;
    std::size_t rows = 0;
    unsigned char fold = 0;
    bool sampled = false;
  };
  std::vector<Spectrum> spectra;
  std::vector<std::size_t> spectrum_of_row(n_rows);
  std::unordered_map<std::uint64_t, std::size_t> index;
  const std::uint64_t salt = mix64(seed ^ 0x6361706669784F46ULL); // "capfixOF"
  for (std::size_t i = 0; i < n_rows; ++i)
  {
    const std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(spec_files[i])) << 32)
                              | static_cast<std::uint32_t>(scans[i]);
    auto it = index.find(key);
    if (it == index.end())
    {
      Spectrum s;
      s.key = key;
      const std::uint64_t h = mix64(salt ^ key);
      s.priority = h;
      s.fold = static_cast<unsigned char>(mix64(h ^ 0xF01DF01DF01DF01DULL) % num_folds);
      it = index.emplace(key, spectra.size()).first;
      spectra.push_back(s);
    }
    spectrum_of_row[i] = it->second;
    ++spectra[it->second].rows;
  }

  std::vector<std::size_t> order(spectra.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&spectra](std::size_t a, std::size_t b)
  {
    if (spectra[a].priority != spectra[b].priority) return spectra[a].priority < spectra[b].priority;
    return spectra[a].key < spectra[b].key;
  });
  OutOfFoldSample out;
  out.num_spectra = spectra.size();
  std::size_t taken = 0;
  for (std::size_t s : order)
  {
    if (spectra[s].rows > max_train - taken) break;
    taken += spectra[s].rows;
    spectra[s].sampled = true;
    ++out.num_sampled_spectra;
  }
  out.fold.resize(n_rows);
  out.sampled_rows.reserve(taken);
  for (std::size_t i = 0; i < n_rows; ++i)
  {
    const Spectrum& s = spectra[spectrum_of_row[i]];
    out.fold[i] = s.fold;
    if (s.sampled) out.sampled_rows.push_back(i);
  }
  return out;
}

}}}  // namespace OpenMS::Internal::Percolator

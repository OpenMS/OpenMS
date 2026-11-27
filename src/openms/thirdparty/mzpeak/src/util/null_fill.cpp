/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <algorithm>
#include <cmath>

#include "mzpeak/util/null_fill.h"

namespace MzPeak::Util {

namespace {

/******************************************************************************/
// Median of an ALREADY-SORTED vector, matching the reference `median`:
// <=2 elements -> the first; odd -> the middle; even -> the mean of the two
// upper-middle elements (mid and mid+1, where mid = n/2).
double median_sorted(const std::vector<double>& d)
{
  std::size_t n = d.size();
  if (n == 0) return 0.0;
  if (n <= 2) return d[0];
  std::size_t mid = n / 2;
  if (n % 2 == 1) return d[mid];
  return (d[mid] + d[mid + 1]) / 2.0;
}

} // namespace

/******************************************************************************/
double predict_delta(const std::vector<double>& beta, double mz)
{
  double acc = 0.0;
  for (std::size_t i = 0; i < beta.size(); ++i) {
    acc += (i == 0) ? beta[i] : beta[i] * std::pow(mz, static_cast<int>(i));
  }
  return acc;
}

/******************************************************************************/
double estimate_median_delta(const std::vector<double>& values,
                             std::size_t begin,
                             std::size_t end)
{
  std::vector<double> deltas;
  if (end > begin + 1) deltas.reserve(end - begin - 1);
  for (std::size_t k = begin + 1; k < end; ++k) {
    deltas.push_back(values[k] - values[k - 1]);
  }
  std::ranges::sort(deltas);
  if (deltas.empty()) return 0.0;

  double m = median_sorted(deltas);
  // Keep only deltas <= the median; the vector stays sorted.
  std::erase_if(deltas, [m](double v) { return v > m; });
  if (deltas.empty()) return m;
  return median_sorted(deltas);
}

/******************************************************************************/
std::vector<double> reconstruct_null_mz(const std::vector<double>& values,
                                        const std::vector<bool>& valid,
                                        const std::vector<double>& beta)
{
  const std::size_t n = values.size();
  // valid is the per-position null bitmap for values; if a caller passes
  // mismatched sizes, indexing valid[i] up to n would read out of bounds.
  // Decline rather than risk UB (callers treat {} as "no reconstruction").
  if (valid.size() != n) return {};

  std::vector<double> out;
  out.reserve(n);

  std::size_t i = 0;
  while (i < n) {
    if (!valid[i]) {
      ++i; // null slots are filled by the adjacent runs
      continue;
    }

    // Maximal run of valid values [s, e).
    std::size_t s = i;
    std::size_t e = s;
    while (e < n && valid[e])
      ++e;

    std::size_t len = e - s;
    double delta = (len > 1) ? estimate_median_delta(values, s, e)
                             : predict_delta(beta, values[s]);

    // Runs are separated by nulls, so a run starting after position 0 has a
    // null immediately before it, and a run ending before n has one after.
    if (s > 0) out.push_back(values[s] - delta);
    for (std::size_t k = s; k < e; ++k)
      out.push_back(values[k]);
    if (e < n) out.push_back(values[e - 1] + delta);

    i = e;
  }

  // Only the interior-paired-null layout reconstructs to the original length.
  // Anything else (e.g. unpaired leading/trailing nulls) is left to the caller.
  if (out.size() != n) return {};
  return out;
}

} // namespace MzPeak::Util

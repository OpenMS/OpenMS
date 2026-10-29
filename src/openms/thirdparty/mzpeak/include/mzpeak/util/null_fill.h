/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <vector>

namespace MzPeak::Util {

/**
 * Evaluate the m/z delta model `beta` at `mz`: a polynomial
 * `beta[0] + beta[1]*mz + beta[2]*mz^2 + ...`.  A single coefficient is a
 * constant model.  Mirrors the reference `MZDeltaModel::predict`.
 */
double predict_delta(const std::vector<double>& beta, double mz);

/**
 * The "median-below-median" delta of the values `[begin, end)`: the median of
 * the consecutive differences that are <= the overall median difference.
 * Mirrors the reference `MedianDeltaEstimator::estimate_median_delta`.
 */
double estimate_median_delta(const std::vector<double>& values,
                             std::size_t begin,
                             std::size_t end);

/**
 * Reconstruct the null-marked m/z values of a profile spectrum.
 *
 * `values[i]` is meaningful only where `valid[i]` is true; the null (invalid)
 * positions are flanking zero-intensity points removed by null-marking.  Each
 * maximal run of valid values emits its values plus one reconstructed value
 * into each adjacent null slot, offset by a local delta (the median spacing of
 * the run, or `predict_delta(beta, value)` for a single-value run).  Mirrors
 * the reference `fill_nulls_for`.
 *
 * Returns the full reconstructed array, or an EMPTY vector if the layout is
 * not the expected interior-paired-null form (the caller should then fall back
 * to leaving nulls as 0).
 */
std::vector<double> reconstruct_null_mz(const std::vector<double>& values,
                                        const std::vector<bool>& valid,
                                        const std::vector<double>& beta);

} // namespace MzPeak::Util

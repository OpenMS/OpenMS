// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "OpenMS/MATH/StatisticFunctions.h"

#include "RunStatistics.h"
#include "Run.h"
#include "Util.h"

// Standard deviation is computed from the interquartile range by
// dividing by 27/20 ≅ 1.349.
const static double IQR_TO_STDDEV = 1.349;

// When a component of the overall score can't be calculated this
// value is used instead to prevent the overall score from going to
// zero.
const static double MIN_SCORE = 3e-7;

namespace OpenMS {
namespace PipEcho {

  /****************************************************************************/
  // NOTE: A lot of this code needs to collect statistics about
  // *identified* features, which we call donors.
  RunStatistics::RunStatistics(const Run& run)
  : log_intensity(init_log_intensity(run)),
    mass_error(init_mass_error(run))
  { }

  /****************************************************************************/
  // Convert a container of numbers into a normal distribution.
  //
  // If `min` is given and the container has fewer elements than
  // `min`, switch to a real standard deviation.  Otherwise convert
  // the interquartile range into a standard deviation.
  //
  // NOTE: `data` will be sorted.
  template <typename T>
  RunStatistics::normal_t RunStatistics::to_normal(T& data,
                                                   const std::optional<double>& min) const
  {
    if (data.size() < 2) return {};

    Math::SummaryStatistics stats(data);
    double stddev = (stats.upperq - stats.lowerq) / IQR_TO_STDDEV;

    if (min.has_value() && data.size() < min) {
      stddev = std::sqrt(stats.variance);
    }

    if (stddev <= 0) return {}; // Guard.
    return boost::math::normal(stats.median, stddev);
  }

  /****************************************************************************/
  RunStatistics::normal_t RunStatistics::init_log_intensity(const Run& run) const {
    std::vector<double> log_intensities(run.donors.storage.size());

    // NOTE: In the PIP-ECHO paper the intensity is divided by the
    // charge.  But in the FlashLFQ implementation that isn't done.
    std::for_each(run.donors.storage.begin(),
                  run.donors.storage.end(),
                  [&](auto& donor) {
                    double intensity = donor->feature.getIntensity();
                    if (intensity > 0) {
                      log_intensities.push_back(std::log2(intensity));
                    }
                  });

    return to_normal(log_intensities, {});
  }

  /****************************************************************************/
  RunStatistics::normal_t RunStatistics::init_mass_error(const Run& run) const {
    std::vector<double> mass_errors(run.donors.storage.size());

    for (auto& donor : run.donors.storage) {
      auto me = Util::feature_mass_error(donor->feature);
      if (me.has_value()) mass_errors.push_back(*me);
    }

    return to_normal(mass_errors, 30);
  }

  /****************************************************************************/
  Score RunStatistics::score(const Feature& donor, const Feature& acceptor) const {
    Score s = {
      .intensity = calc_intensity_score(acceptor),
      .rt_diff_error = acceptor.getRT() - donor.getRT(),
      .mass_error = calc_mass_error_score(acceptor),
      .mbr_score = MIN_SCORE
    };

    // IMPORTANT: This is the number of individual scores that are
    // multiplied below.  Please update this if you change the
    // following call to `std::pow`.
    const static double measures = 3.0;

    s.mbr_score =
      100.0 *
      std::pow(s.intensity * s.rt_diff_error * s.mass_error,
               1.0/measures);

    return s;
  }

  /****************************************************************************/
  double RunStatistics::calc_intensity_score(const Feature& acceptor) const {
    return calc_score_using(log_intensity,
                            std::log2(acceptor.getIntensity()));
  }

  /****************************************************************************/
  double RunStatistics::calc_mass_error_score(const Feature& acceptor) const {
    auto me = Util::feature_mass_error(acceptor);
    if (!me.has_value()) return MIN_SCORE;
    return calc_score_using(mass_error, *me);
  }

  /****************************************************************************/
  double RunStatistics::calc_score_using(const normal_t& dist, double value) const {
    if (!dist.has_value()) return MIN_SCORE;
    double diff = std::fabs(dist->mean() - value);
    return 2 * boost::math::cdf(*dist, dist->mean() - diff);
  }
}}

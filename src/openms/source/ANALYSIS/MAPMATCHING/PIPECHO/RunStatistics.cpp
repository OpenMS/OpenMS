// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include "Run.h"
#include "RunStatistics.h"
#include "Util.h"

#include <algorithm>
#include <cmath>

// Standard deviation is computed from the interquartile range by
// dividing by 27/20 ≅ 1.349.
const static double IQR_TO_STDDEV = 1.349;

// When a component of the overall score can't be calculated this
// value is used instead to prevent the overall score from going to
// zero.
const static double MIN_SCORE = 3e-7;

namespace OpenMS::PipEcho
{

// File-local helper (internal linkage); not part of the module's header API.
namespace
{
/******************************************************************************/
/// Ensure the given value is within bounds (i.e. not NaN.)
double bound(double d)
{
  if (std::isnan(d)) { return MIN_SCORE; }
  else
  {
    return d;
  }
}
} // namespace

/******************************************************************************/
// NOTE: A lot of this code needs to collect statistics about
// *identified* features, which we call donors.
RunStatistics::RunStatistics(const Run& run, bool enable_im):
    log_intensity(init_log_intensity(run)),
    mass_error(init_mass_error(run)),
    im_tolerance(init_im_tolerance(run, enable_im))
{
}

/******************************************************************************/
// Convert a container of numbers into a normal distribution.
//
// If `min` is given and the container has fewer elements than
// `min`, switch to a real standard deviation.  Otherwise convert
// the interquartile range into a standard deviation.
//
// NOTE: `data` will be sorted.
template<typename T>
RunStatistics::normal_t
RunStatistics::to_normal(T& data, const std::optional<double>& min) const
{
  if (data.size() < 2) return {};

  Math::SummaryStatistics stats(data);
  double stddev = (stats.upperq - stats.lowerq) / IQR_TO_STDDEV;

  if (min.has_value() && data.size() < min)
  {
    stddev = std::sqrt(stats.variance);
  }

  if (stddev <= 0) return {}; // Guard.
  return boost::math::normal(stats.median, stddev);
}

/******************************************************************************/
RunStatistics::normal_t RunStatistics::init_log_intensity(const Run& run) const
{
  std::vector<double> log_intensities;
  log_intensities.reserve(run.donors.storage.size());

  // NOTE: In the PIP-ECHO paper the intensity is divided by the
  // charge.  But in the FlashLFQ implementation that isn't done.
  std::for_each(run.donors.storage.begin(), run.donors.storage.end(),
                [&](auto& donor) {
                  double intensity = donor->feature.getIntensity();
                  if (intensity > 0)
                  {
                    log_intensities.push_back(std::log2(intensity));
                  }
                });

  return to_normal(log_intensities, {});
}

/******************************************************************************/
RunStatistics::normal_t RunStatistics::init_mass_error(const Run& run) const
{
  std::vector<double> mass_errors;
  mass_errors.reserve(run.donors.storage.size());

  for (auto& donor : run.donors.storage)
  {
    auto me = Util::feature_mass_error(donor->feature);
    if (me.has_value()) mass_errors.push_back(*me);
  }

  return to_normal(mass_errors, 30);
}

/******************************************************************************/
RunStatistics::normal_t RunStatistics::init_im_tolerance(const Run& run, bool enable_im) const
{
  // Ion mobility is enabled globally (all runs) or not at all; when disabled,
  // build no tolerance so the IM feature is omitted from every run's scores.
  if (! enable_im) return {};

  // Derive an ion-mobility tolerance from the data, using the IM elution
  // width (IM_max - IM_min) of identified features as a proxy for the scale on
  // which the same ion's mobility may vary. Half the typical width is used as
  // the standard deviation of a zero-mean tolerance distribution. When the run
  // carries no ion mobility (no widths), this returns nullopt and the IM
  // feature is silently dropped.
  std::vector<double> widths;
  widths.reserve(run.donors.storage.size());

  for (auto& donor : run.donors.storage)
  {
    auto w = Util::feature_im_width(donor->feature);
    if (w.has_value()) widths.push_back(*w);
  }

  if (widths.size() < 2) return {};

  // Use the median width directly (Math::median is well-defined for >= 2
  // values, unlike Math::SummaryStatistics whose quantiles hit an empty
  // sub-range at exactly 2 elements).
  double stddev = Math::median(widths.begin(), widths.end(), false) / 2.0;

  if (stddev <= 0) return {}; // Guard.
  return boost::math::normal(0.0, stddev);
}

/******************************************************************************/
Score RunStatistics::score(const Feature& donor, const Feature& acceptor,
                           std::optional<double> donor_rt_override,
                           RtScoreMode rt_mode,
                           const boost::math::normal* rt_dist) const
{
  // Neutral RT agreement used when a calibrated mode is requested but the run
  // pair has no estimable RT prediction-error distribution (too few anchors).
  // Codex review: an absent distribution must NOT become a hard penalty, or it
  // would inject a spurious "this pair is unreliable" axis that the SVM can pick
  // up (and that the decoy null does not symmetrically reproduce); a neutral
  // value keeps the column rectangular without confounding the FDR.
  constexpr double NEUTRAL_RT = 0.5;

  // For decoys the RT difference is measured against the randomised RT anchor,
  // not the donor's true RT (otherwise decoys carry an artificially huge RT
  // error and the classifier separates them trivially -- a label leak). On the
  // local adaptive RT path donor_rt is the locally PREDICTED acceptor RT centre,
  // so the SIGNED residual below is exactly what the RT prediction-error
  // distribution (built leave-one-out from anchors) is calibrated against.
  const double donor_rt = donor_rt_override.value_or(donor.getRT());
  const double signed_rt = acceptor.getRT() - donor_rt;
  Score s = {.intensity = calc_intensity_score(acceptor),
             .rt_diff_error = std::fabs(signed_rt),
             .mass_error = calc_mass_error_score(donor, acceptor),
             .im_diff_score = -1.0, // sentinel: ion mobility not applicable
             .isotope_score = calc_isotope_score(donor, acceptor),
             .rt_score = -1.0,      // sentinel: RT not calibrated (Raw path)
             .mbr_score = MIN_SCORE};

  // Calibrate the RT residual into a [0, 1] agreement score (FlashLFQ-style:
  // the two-tailed CDF of the signed residual under the run pair's RT
  // prediction-error distribution). Only on the local adaptive RT path; the
  // global/baseline path keeps the raw rt_diff_error and the -1 sentinel, so its
  // output stays byte-identical.
  if (rt_mode != RtScoreMode::Raw)
  {
    s.rt_score = (rt_dist != nullptr) ? calc_score_using(*rt_dist, signed_rt)
                                      : NEUTRAL_RT;
  }

  // The MBR bootstrap score is the geometric mean of the agreement scores.
  // Retention time was historically EXCLUDED because the raw rt_diff_error is an
  // un-normalised |Δrt| with inverted monotonicity (a larger RT error would
  // *increase* the score). With a CALIBRATED rt_score on the local path that
  // objection is gone, so RtScoreMode::SvmAndMbr folds it in like the other
  // agreement scores. Ion mobility is an optional extra factor, included only
  // when the data carries it (e.g. timsTOF/.d). `measures` tracks how many scores
  // are multiplied below; keep it in sync with the product.
  double product = s.intensity * s.mass_error;
  double measures = 2.0;

  if (auto im = calc_im_score(donor, acceptor); im.has_value())
  {
    s.im_diff_score = *im;
    product *= *im;
    measures += 1.0;
  }

  // Fold the calibrated RT score into the geometric mean only when it was
  // actually calibrated against a real distribution. For a sparse pair with no
  // RT distribution rt_score is the NEUTRAL 0.5; multiplying that in (and bumping
  // the exponent) would NOT be neutral for mbr_score, so skip it there -- the
  // 0.5 still goes to the SVM column for rectangularity (Codex review).
  if (rt_mode == RtScoreMode::SvmAndMbr && rt_dist != nullptr)
  {
    product *= s.rt_score;
    measures += 1.0;
  }

  s.mbr_score = bound(100.0 * std::pow(product, 1.0 / measures));

  return s;
}

/******************************************************************************/
double RunStatistics::calc_intensity_score(const Feature& acceptor) const
{
  return calc_score_using(log_intensity, std::log2(acceptor.getIntensity()));
}

/******************************************************************************/
double RunStatistics::calc_mass_error_score(const Feature& donor,
                                            const Feature& acceptor) const
{
  // The acceptor is unidentified by construction, so its own ID cannot give a
  // mass error.  Compute the acceptor m/z error against the donor's peptide
  // sequence/charge -- the identity that would be transferred to it.
  auto donor_hit = Util::feature_hit(donor);
  if (! donor_hit.has_value()) return MIN_SCORE;
  double theoretical = donor_hit->getSequence().getMZ(donor_hit->getCharge());
  return calc_score_using(mass_error, Math::getPPM(acceptor.getMZ(), theoretical));
}

/******************************************************************************/
std::optional<double>
RunStatistics::calc_im_score(const Feature& donor, const Feature& acceptor) const
{
  // Not applicable: ion mobility is globally disabled (no tolerance was built),
  // so the IM feature is omitted from the geometric mean for EVERY pair.
  if (! im_tolerance.has_value()) return std::nullopt;

  auto donor_im = Util::feature_ion_mobility(donor);
  auto acceptor_im = Util::feature_ion_mobility(acceptor);

  // IM is globally enabled but one of these features lacks an IM annotation
  // (anomalous in IM data). Penalise with MIN_SCORE rather than dropping the
  // feature, so `measures` (and the SVM IM column) stay consistent across all
  // pairs -- a per-pair drop would put this acceptor's mbr_score on a different
  // scale than its peers and desync the global use_im decision in TransferFDRModel.
  if (! donor_im.has_value() || ! acceptor_im.has_value()) return MIN_SCORE;

  // Score the IM difference under the zero-mean tolerance: a small difference
  // (same ion) scores near 1, a large one decays towards MIN_SCORE. Reuses the
  // two-tailed CDF helper (mean 0 => 2 * cdf(N(0, sigma), -|delta|)).
  return calc_score_using(im_tolerance, *acceptor_im - *donor_im);
}

/******************************************************************************/
double RunStatistics::calc_isotope_score(const Feature& donor, const Feature& acceptor) const
{
  // Neutral "uncorrelated" value when the envelope is unavailable/too short.
  constexpr double NA = 0.5;

  auto obs = Util::feature_obs_envelope(acceptor);
  auto donor_hit = Util::feature_hit(donor);
  if (! obs.has_value() || ! donor_hit.has_value()) { return NA; }

  // Donor peptide's theoretical isotope envelope (neutral formula -- the relative
  // isotope intensities are charge-independent), to the observed envelope length.
  IsotopeDistribution iso = donor_hit->getSequence().getFormula().getIsotopeDistribution(
    CoarseIsotopePatternGenerator(obs->size()));

  std::vector<double> theo;
  theo.reserve(iso.size());
  for (const auto& peak : iso) { theo.push_back(peak.getIntensity()); }

  const std::size_t n = std::min(obs->size(), theo.size());
  if (n < 3) { return NA; } // Pearson is degenerate for < 3 points

  const double r = Math::pearsonCorrelationCoefficient(
    obs->begin(), obs->begin() + n, theo.begin(), theo.begin() + n);
  if (! std::isfinite(r)) { return NA; } // a constant envelope yields NaN

  return (r + 1.0) / 2.0; // map [-1, 1] -> [0, 1] (1 = identical envelope shape)
}

/******************************************************************************/
double RunStatistics::calc_score_using(const normal_t& dist, double value) const
{
  if (! dist.has_value()) return MIN_SCORE;
  return calc_score_using(*dist, value);
}

/******************************************************************************/
double RunStatistics::calc_score_using(const boost::math::normal& dist,
                                       double value) const
{
  double diff = std::fabs(dist.mean() - value);
  return bound(2 * boost::math::cdf(dist, dist.mean() - diff));
}

} // namespace OpenMS::PipEcho

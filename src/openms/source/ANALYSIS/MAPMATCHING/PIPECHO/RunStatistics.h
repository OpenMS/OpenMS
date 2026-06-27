// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "Run.h"
#include "Score.h"

#include <boost/math/distributions/normal.hpp>

namespace OpenMS::PipEcho
{

/******************************************************************************/
class RunStatistics
{
public:
  /// Generate statistics for the given Run object.
  ///
  /// @param enable_im Whether ion mobility is used as a scoring feature. This
  /// is a global, all-or-nothing decision made by the caller across all runs
  /// (so the geometric-mean MBR score stays on one scale); when false, no IM
  /// tolerance is built and the IM feature is omitted everywhere.
  RunStatistics(const Run&, bool enable_im);

  /// Score the given donor and acceptor.
  ///
  /// For decoy candidates the donor is hypothetically eluting at a randomised
  /// (wrong) retention time; pass it via @p donor_rt_override so the
  /// retention-time-difference feature is measured against that randomised RT
  /// rather than the donor's true RT (PIP-ECHO: the RT feature is the deviation
  /// from the predicted/randomised RT). Only the RT feature is affected; mass
  /// error and intensity remain properties of the donor identity / acceptor.
  Score score(const Feature&, const Feature&,
              std::optional<double> donor_rt_override = std::nullopt) const;

private:
  // Normal distribution: <mean, stddev>.
  using normal_t = std::optional<boost::math::normal>;

  template<typename T>
  normal_t to_normal(T&, const std::optional<double>&) const;

  normal_t init_log_intensity(const Run&) const;
  normal_t init_mass_error(const Run&) const;
  normal_t init_im_tolerance(const Run&, bool enable_im) const;
  double calc_score_using(const normal_t&, double) const;
  double calc_intensity_score(const Feature&) const;
  double calc_mass_error_score(const Feature& donor, const Feature& acceptor) const;

  /// Ion-mobility agreement score for a donor/acceptor pair.
  ///
  /// Returns nullopt when ion mobility is not applicable -- the run carries no
  /// ion-mobility data, or either feature lacks an IM annotation -- so the
  /// caller can omit the feature entirely (rather than penalise the score).
  /// Otherwise returns a value in [0, 1] (1 = identical mobility). For decoys
  /// the donor's real IM is used (IM is independent of the randomised decoy
  /// RT), mirroring the mass-error feature.
  std::optional<double> calc_im_score(const Feature& donor, const Feature& acceptor) const;

  /// Isotope-distribution agreement score for a donor/acceptor pair: the Pearson
  /// correlation between the acceptor's observed isotope envelope and the donor
  /// peptide's theoretical envelope, mapped to [0, 1]. Returns 0.5
  /// ("uncorrelated") when the envelope is unavailable or too short. Identity-
  /// coupled (compares the acceptor signal to the DONOR's theoretical pattern),
  /// so a coincidental m/z match scores low.
  double calc_isotope_score(const Feature& donor, const Feature& acceptor) const;

private:
  // Log intensity distribution.
  normal_t log_intensity;

  // Mass error (PPM) distribution.
  normal_t mass_error;

  // Ion-mobility tolerance, modelled as a half-normal centred at 0 (a zero-mean
  // normal scored two-tailed). Present only when the run carries ion mobility.
  normal_t im_tolerance;
};

} // namespace OpenMS::PipEcho

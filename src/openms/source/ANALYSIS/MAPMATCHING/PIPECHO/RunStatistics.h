// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
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
  RunStatistics(const Run&);

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
  double calc_score_using(const normal_t&, double) const;
  double calc_intensity_score(const Feature&) const;
  double calc_mass_error_score(const Feature& donor, const Feature& acceptor) const;

private:
  // Log intensity distribution.
  normal_t log_intensity;

  // Mass error (PPM) distribution.
  normal_t mass_error;
};

} // namespace OpenMS::PipEcho

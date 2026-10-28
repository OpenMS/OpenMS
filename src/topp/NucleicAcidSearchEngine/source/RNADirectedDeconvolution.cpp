// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Sam Wein $
// $Authors: Sam Wein $
// --------------------------------------------------------------------------
#include "RNADirectedDeconvolution.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <set>

namespace OpenMS
{
  void RNADirectedDeconvolution::deconvolve(const MSSpectrum& spectrum, int scan_number,
    const std::vector<std::pair<double, int>>& queries)
  {
    queries_.clear();
    if (!queries.empty() && spectrum.getMSLevel() != 1)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor queries require an MS1 spectrum.");
    }
    for (const auto& query : queries)
    {
      if (!std::isfinite(query.first) || query.first <= 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Query m/z must be finite and positive.");
      }
    }
    try
    {
      for (const auto& query : queries)
      {
        Precursor precursor;
        precursor.setMZ(query.first);
        precursor.setCharge(query.second);
        queries_.push_back(precursor);
      }
      performSpectrumDeconvolution(spectrum, scan_number, PeakGroup());
    }
    catch (...)
    {
      queries_.clear();
      throw;
    }
    queries_.clear();
  }

  std::optional<SpectralDeconvolution::MassRanges> RNADirectedDeconvolution::selectCandidateMassRanges_(
    const CandidateSelectionContext& context, const DeconvolvedSpectrum* evidence_groups,
    const std::vector<Size>& valid) const
  {
    if (queries_.empty()) return std::nullopt;
    const double carrier = FLASHHelperClasses::getChargeMass(context.positive);
    // Cover initial isotope recruitment, apex reassignment, and all ten offset
    // refinement steps. This bound is calculated from the configured averagine,
    // rather than truncating the precursor search to a fixed isotope list.
    const auto isotope_span = [&](double mass)
    {
      mass = std::min(context.max_mass, mass + context.averagine.getMaxIsotopeIndex() * context.isotope_distance);
      return int(context.averagine.getLeftCountFromApex(mass) + context.averagine.getRightCountFromApex(mass) +
        context.averagine.getApexIndex(mass) + 10 * (context.averagine.getApexIndex(mass) / 4 + 1) + 4);
    };
    std::vector<std::pair<double, double>> intervals;
    for (const auto& precursor : queries_)
    {
      if (!std::isfinite(precursor.getMZ()) || precursor.getMZ() <= 0) continue;
      const int reported = std::abs(precursor.getCharge());
      for (int z = context.min_charge; z <= context.max_charge; ++z)
      {
        if (reported != 0 && reported != z) continue;
        const double mass = (precursor.getMZ() - carrier) * z;
        // Initial apex bins can shift continuously when overlapping charge
        // states refine the mass. Restricting them to exact isotope-spaced
        // bands can discard a competing envelope before global filtering.
        const double delta = isotope_span(mass) * context.isotope_distance + 4 * context.tolerance * mass;
        intervals.emplace_back(mass - context.averagine.getLastIndex(mass) * context.isotope_distance - delta, mass + delta);
      }
    }
    if (!evidence_groups) return intervals;
    // Charge-error rejection compares groups that share an observed signal peak.
    // Include every alternative charge interpretation of the queried groups'
    // evidence before applying the ordinary global filters. Do not mark these
    // groups as targeted: they receive no exemption from competition.
    std::set<double> evidence;
    for (Size i : valid)
    {
      const auto& group = (*evidence_groups)[i];
      bool relevant = false;
      for (const auto& precursor : queries_)
      {
        const double delta = precursor.getMZ() * context.tolerance;
        for (const auto& peak : group)
        {
          if ((precursor.getCharge() == 0 || std::abs(precursor.getCharge()) == peak.abs_charge) &&
              peak.isotopeIndex >= 0 && std::abs(peak.mz - precursor.getMZ()) <= delta)
          {
            relevant = true;
            break;
          }
        }
        if (relevant) break;
      }
      if (relevant)
      {
        for (const auto& peak : group) evidence.insert(peak.mz);
      }
    }
    // Also retain the seed intervals for neighboring groups and isotope aliases.
    // Merge overlapping profile-sample bands before expanding charge/isotope
    // aliases; otherwise thousands of samples create millions of duplicate bands.
    std::vector<std::pair<double, double>> observed_ranges;
    for (double mz : evidence)
    {
      if (observed_ranges.empty() || mz - observed_ranges.back().second >
          2 * context.tolerance * mz)
      {
        observed_ranges.emplace_back(mz, mz);
      }
      else observed_ranges.back().second = mz;
    }
    for (const auto& [low_mz, high_mz] : observed_ranges)
    {
      for (int z = 1; z <= context.max_charge; ++z)
      {
        const double mass = ((low_mz + high_mz) / 2 - carrier) * z;
        const double delta = isotope_span(mass) * context.isotope_distance +
          4 * context.tolerance * mass + (high_mz - low_mz) * z / 2;
        intervals.emplace_back(mass - delta, mass + delta);
      }
    }
    return intervals;
  }
}

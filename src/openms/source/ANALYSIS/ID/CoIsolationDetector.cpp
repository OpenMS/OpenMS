// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/CoIsolationDetector.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace OpenMS
{
  const std::string CoIsolationDetector::META_MZ = "coisolated_mz";
  const std::string CoIsolationDetector::META_CHARGE = "coisolated_charge";

  namespace
  {
    constexpr double NO_ENVELOPE_KL = 10.0;
    constexpr Size MAX_MONO_WALK = 3;
    constexpr int MAX_ISOTOPE_STEPS = 6;

    /// most intense peak within [lo, hi] of a sorted spectrum, 0 if none
    double maxIntensityIn(const MSSpectrum& s, double lo, double hi)
    {
      double best = 0.0;
      for (auto it = s.MZBegin(lo); it != s.end() && it->getMZ() <= hi; ++it)
      {
        best = std::max(best, static_cast<double>(it->getIntensity()));
      }
      return best;
    }

    bool hasPeakIn(const MSSpectrum& s, double lo, double hi)
    {
      auto it = s.MZBegin(lo);
      return it != s.end() && it->getMZ() <= hi;
    }
  }

  double CoIsolationDetector::isotopeKL(const MSSpectrum& ms1, double mono_mz, int charge, double tol_da, Size n_isotopes)
  {
    if (ms1.empty() || n_isotopes == 0 || charge <= 0) return NO_ENVELOPE_KL;
    const double spacing = Constants::C13C12_MASSDIFF_U / charge;
    std::vector<double> observed(n_isotopes, 0.0);
    for (Size k = 0; k < n_isotopes; ++k)
    {
      const double target = mono_mz + static_cast<double>(k) * spacing;
      observed[k] = maxIntensityIn(ms1, target - tol_da, target + tol_da);
    }
    const double sum = std::accumulate(observed.begin(), observed.end(), 0.0);
    if (sum <= 0.0) return NO_ENVELOPE_KL;

    const double neutral_mass = (mono_mz - Constants::PROTON_MASS_U) * charge;
    CoarseIsotopePatternGenerator generator(n_isotopes);
    const IsotopeDistribution averagine = generator.estimateFromPeptideWeight(neutral_mass);
    std::vector<double> theo(n_isotopes, 0.0);
    for (Size k = 0; k < n_isotopes && k < averagine.size(); ++k) theo[k] = averagine[k].getIntensity();
    const double theo_sum = std::accumulate(theo.begin(), theo.end(), 0.0);

    constexpr double eps = 1e-9;
    double kl = 0.0;
    for (Size k = 0; k < n_isotopes; ++k)
    {
      const double o = observed[k] / sum;
      if (o <= 0.0) continue;
      const double t = theo_sum > 0.0 ? theo[k] / theo_sum : 0.0;
      kl += o * std::log(o / std::max(t, eps));
    }
    return kl;
  }

  std::vector<CoIsolationDetector::CoIsolated> CoIsolationDetector::detect(const MSSpectrum& ms1,
                                                                          double window_lower,
                                                                          double window_upper,
                                                                          double selected_mz,
                                                                          int min_charge,
                                                                          int max_charge,
                                                                          double tol_da,
                                                                          double max_kl,
                                                                          Size max_n)
  {
    std::vector<CoIsolated> out;
    if (max_n == 0 || ms1.empty() || window_upper <= window_lower) return out;
    min_charge = std::max(min_charge, 1);

    // peaks in the window, most intense first
    std::vector<std::pair<double, double>> candidates; // (m/z, intensity)
    for (auto it = ms1.MZBegin(window_lower); it != ms1.end() && it->getMZ() <= window_upper; ++it)
    {
      candidates.emplace_back(it->getMZ(), it->getIntensity());
    }
    std::stable_sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return a.second > b.second; });

    // a peak of the selected precursor's envelope (either side: the instrument may have selected M+1)
    auto is_selected = [&](double mz)
    {
      const double d = std::fabs(mz - selected_mz);
      if (d <= tol_da) return true;
      for (int z = min_charge; z <= max_charge; ++z)
      {
        for (int k = 1; k < MAX_ISOTOPE_STEPS; ++k)
        {
          if (std::fabs(d - k * Constants::C13C12_MASSDIFF_U / z) <= tol_da) return true;
        }
      }
      return false;
    };

    for (const auto& candidate : candidates)
    {
      // not a structured binding: clang cannot capture one in a lambda when OpenMP is enabled
      const double mz = candidate.first;
      if (out.size() >= max_n) break;
      if (is_selected(mz)) continue;
      // an isotope peak of an envelope already accepted
      const bool known = std::any_of(out.begin(), out.end(), [&](const CoIsolated& c)
      {
        const double d = std::fabs(mz - c.mono_mz);
        for (int k = 0; k < MAX_ISOTOPE_STEPS; ++k)
        {
          if (std::fabs(d - k * Constants::C13C12_MASSDIFF_U / c.charge) <= tol_da) return true;
        }
        return false;
      });
      if (known) continue;

      // the charge whose envelope fits averagine best
      CoIsolated best;
      best.kl = std::numeric_limits<double>::max();
      for (int z = min_charge; z <= max_charge; ++z)
      {
        const double kl = isotopeKL(ms1, mz, z, tol_da);
        if (kl <= max_kl && kl < best.kl) best = {mz, z, kl};
      }
      if (best.charge == 0) continue;

      // move down to the monoisotopic peak
      for (Size step = 0; step < MAX_MONO_WALK; ++step)
      {
        const double lower = best.mono_mz - Constants::C13C12_MASSDIFF_U / best.charge;
        if (!hasPeakIn(ms1, lower - tol_da, lower + tol_da) || is_selected(lower)) break;
        const double lower_kl = isotopeKL(ms1, lower, best.charge, tol_da);
        if (lower_kl >= best.kl) break;
        best.mono_mz = lower;
        best.kl = lower_kl;
      }
      out.push_back(best);
    }
    return out;
  }

  Size CoIsolationDetector::annotate(PeakMap& exp, int min_charge, int max_charge, double tolerance, bool tolerance_ppm,
                                     double max_kl, Size max_n)
  {
    Size found = 0;
    for (MSSpectrum& s : exp)
    {
      if (s.getMSLevel() == 1 && !s.isSorted()) s.sortByPosition();
    }
    const PeakMap& const_exp = exp;
    for (auto it = exp.begin(); it != exp.end(); ++it)
    {
      if (it->getMSLevel() != 2 || it->getPrecursors().empty()) continue;
      const Precursor& precursor = it->getPrecursors()[0];
      const double lower_offset = precursor.getIsolationWindowLowerOffset();
      const double upper_offset = precursor.getIsolationWindowUpperOffset();
      if (lower_offset <= 0.0 && upper_offset <= 0.0) continue;
      const auto survey = const_exp.getPrecursorSpectrum(const_exp.begin() + (it - exp.begin()));
      if (survey == const_exp.end() || survey->getMSLevel() != 1) continue;

      double center = precursor.getMZ();
      if (precursor.metaValueExists("isolation window target m/z"))
      {
        const double target = precursor.getMetaValue("isolation window target m/z");
        if (target > 0.0) center = target;
      }
      const MSSpectrum& ms1 = *survey;
      const double tol_da = std::max(tolerance_ppm ? precursor.getMZ() * tolerance * 1e-6 : tolerance, MIN_TOLERANCE_DA);
      const auto coisolated = detect(ms1, center - lower_offset, center + upper_offset, precursor.getMZ(),
                                     min_charge, max_charge, tol_da, max_kl, max_n);
      DoubleList mzs;
      IntList charges;
      for (const auto& c : coisolated)
      {
        mzs.push_back(c.mono_mz);
        charges.push_back(c.charge);
      }
      it->setMetaValue(META_MZ, mzs);
      it->setMetaValue(META_CHARGE, charges);
      found += coisolated.size();
    }
    return found;
  }
} // namespace OpenMS

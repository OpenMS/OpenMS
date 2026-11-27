// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "NuXLIonLikelihood.h"

#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>
#include <cmath>

namespace OpenMS
{
  NuXLIonLikelihood::NuXLIonLikelihood(double tolerance, bool ppm, int max_fragment_charge) :
    tolerance_(tolerance),
    ppm_(ppm),
    max_fragment_charge_(max_fragment_charge)
  {
    Param p = tsg_.getParameters();
    p.setValue("add_b_ions", "true");
    p.setValue("add_y_ions", "true");
    p.setValue("add_a_ions", "false");
    p.setValue("add_first_prefix_ion", "true");
    p.setValue("add_losses", "false");
    p.setValue("isotope_model", "none");
    p.setValue("add_precursor_peaks", "false");
    p.setValue("add_abundant_immonium_ions", "false");
    p.setValue("add_metainfo", "true");
    tsg_.setParameters(p);
  }

  void NuXLIonLikelihood::setSpectra(const PeakMap& spectra)
  {
    peaks_.reset(spectra.size());
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i)
    {
      peaks_.assign(static_cast<Size>(i), spectra[i]);
    }
  }

  AASequence NuXLIonLikelihood::reversedSequence(const AASequence& peptide)
  {
    const Size n = peptide.size();
    if (n < 3) return peptide;
    AASequence reversed;
    for (Size i = n - 1; i-- > 0;) { reversed += &peptide[i]; }
    reversed += &peptide[n - 1];
    if (peptide.hasNTerminalModification()) { reversed.setNTerminalModification(peptide.getNTerminalModification()); }
    if (peptide.hasCTerminalModification()) { reversed.setCTerminalModification(peptide.getCTerminalModification()); }
    return reversed;
  }

  NuXLIonLikelihood::Outcomes NuXLIonLikelihood::matchHypothesis_(Size scan, const AASequence& peptide, int precursor_charge,
                                                                  const std::vector<double>& shifts, int site) const
  {
    Outcomes out;
    if (scan >= peaks_.size() || peaks_.peaks(scan) == 0 || peptide.size() < 2) return out;
    const int max_charge = std::max(1, std::min(max_fragment_charge_, precursor_charge - 1));
    PeakSpectrum theo;
    tsg_.getSpectrum(theo, peptide, 1, max_charge);
    const Model& matcher = models_[0][0]; // matching only depends on the context set, which all models share

    if (shifts.empty())
    {
      matcher.matchIons(peaks_, scan, theo, peptide, precursor_charge, tolerance_, ppm_, out.ions);
      out.groups.assign(out.ions.size(), LINEAR);
      return out;
    }

    // Ions that carry the fragment adduct: those containing the cross-linked residue, or all if it is not localized.
    // matchIons reports one ion per parseable theoretical peak in spectrum order, so the k-th ion of every
    // alternative below is the same theoretical ion.
    const Size n = peptide.size();
    const auto& names = theo.getStringDataArrays()[0];
    const auto& charges = theo.getIntegerDataArrays()[0];
    std::vector<Size> parsed;
    std::vector<char> carries_adduct;
    for (Size i = 0; i < theo.size(); ++i)
    {
      bool prefix = false;
      Size ordinal = 0;
      if (!Model::parseIonName(names[i], prefix, ordinal) || ordinal >= n) continue;
      parsed.push_back(i);
      const SignedSize s = site;
      carries_adduct.push_back(site < 0 ? 1 : (prefix ? s < static_cast<SignedSize>(ordinal) : s >= static_cast<SignedSize>(n - ordinal)));
    }

    // alternatives: every fragment adduct, and for an unlocalized cross-link also no adduct
    std::vector<Model::Ion> alternative;
    bool first = true;
    auto merge = [&](const PeakSpectrum& spectrum)
    {
      matcher.matchIons(peaks_, scan, spectrum, peptide, precursor_charge, tolerance_, ppm_, alternative);
      if (first)
      {
        out.ions = alternative;
        first = false;
        return;
      }
      // best outcome (lowest rank bin, then lowest error bin; absence is the last outcome) over the alternatives;
      // the complementary ion counts as matched if it matched in any alternative (lowest bit of the RICH context)
      static_assert(Model::COMPLEMENT_STATES == 2);
      for (Size k = 0; k < out.ions.size(); ++k)
      {
        out.ions[k].outcome = std::min(out.ions[k].outcome, alternative[k].outcome);
        out.ions[k].context |= (alternative[k].context & 1u);
      }
    };

    if (site < 0) merge(theo);
    PeakSpectrum shifted = theo;
    for (double shift : shifts)
    {
      for (Size k = 0; k < parsed.size(); ++k)
      {
        const Size i = parsed[k];
        shifted[i].setMZ(carries_adduct[k] ? theo[i].getMZ() + shift / static_cast<double>(std::max(1, charges[i])) : theo[i].getMZ());
      }
      merge(shifted);
    }

    out.groups.resize(out.ions.size());
    for (Size k = 0; k < out.ions.size(); ++k)
    {
      out.groups[k] = site < 0 ? XL_UNLOCALIZED : (carries_adduct[k] ? XL_SHIFTED : XL_UNSHIFTED);
    }
    return out;
  }

  NuXLIonLikelihood::Evidence NuXLIonLikelihood::match(Size scan, const AASequence& peptide, int precursor_charge,
                                                       const std::vector<double>& shifts, int site) const
  {
    Evidence e;
    e.signal = matchHypothesis_(scan, peptide, precursor_charge, shifts, site);
    const AASequence reversed = reversedSequence(peptide);
    if (reversed == peptide) return e;
    const int n = static_cast<int>(peptide.size());
    const int reversed_site = (site >= 0 && site < n - 1) ? n - 2 - site : site; // the cross-linked residue moves with the reversal
    e.noise = matchHypothesis_(scan, reversed, precursor_charge, shifts, reversed_site);
    return e;
  }

  void NuXLIonLikelihood::resetFold(Size fold)
  {
    for (Model& m : models_[fold]) { m = Model(); }
  }

  void NuXLIonLikelihood::addTraining(Size fold, const Evidence& evidence)
  {
    std::array<std::vector<Model::Ion>, GROUPS> by_group;
    auto add = [&](const Outcomes& o, bool noise)
    {
      for (auto& g : by_group) g.clear();
      for (Size k = 0; k < o.ions.size(); ++k) by_group[o.groups[k]].push_back(o.ions[k]);
      for (Size g = 0; g < GROUPS; ++g)
      {
        if (!by_group[g].empty()) models_[fold][g].addObservations(by_group[g], noise);
      }
    };
    add(evidence.signal, false);
    add(evidence.noise, true);
  }

  void NuXLIonLikelihood::finalize()
  {
    for (auto& fold : models_)
    {
      for (Model& m : fold) { m.finalize(); }
    }
  }

  NuXLIonLikelihood::Features NuXLIonLikelihood::score(Size fold, const Outcomes& outcomes) const
  {
    Features f;
    if (outcomes.ions.empty()) return f;
    const Size absent = models_[fold][0].absentOutcome();
    double presence_total = 0.0, presence_matched = 0.0;
    std::vector<std::pair<double, bool>> predicted;
    predicted.reserve(outcomes.ions.size());
    for (Size k = 0; k < outcomes.ions.size(); ++k)
    {
      const Model& m = models_[fold][outcomes.groups[k]];
      const Model::Ion& ion = outcomes.ions[k];
      const double llr = m.logLikelihoodRatio(ion.context, ion.outcome);
      f.llr += llr;
      if (outcomes.groups[k] == XL_SHIFTED || outcomes.groups[k] == XL_UNLOCALIZED) f.llr_shifted += llr;
      const double presence = m.presenceProbability(ion.context);
      const bool matched = ion.outcome != absent;
      presence_total += presence;
      if (matched) presence_matched += presence;
      predicted.emplace_back(presence, matched);
    }
    f.explained_presence = presence_total > 0.0 ? presence_matched / presence_total : 0.0;
    std::stable_sort(predicted.begin(), predicted.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    const Size top_k = std::min<Size>(6, predicted.size());
    const Size observed = static_cast<Size>(std::count_if(predicted.begin(), predicted.begin() + top_k, [](const auto& p) { return p.second; }));
    f.top_predicted_observed = static_cast<double>(observed) / static_cast<double>(top_k);
    return f;
  }
} // namespace OpenMS

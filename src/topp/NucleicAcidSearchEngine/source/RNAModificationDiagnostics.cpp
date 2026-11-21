// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "RNAModificationDiagnostics.h"

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>

namespace OpenMS
{
bool RNAModificationDiagnostics::Evidence::accepts(const NASequence& candidate) const
{
  return std::all_of(candidate.begin(), candidate.end(),
                     [&](const auto& residue) { return ! tested.contains(residue.getCode()) || allowed.contains(residue.getCode()); });
}

RNAModificationDiagnostics::RNAModificationDiagnostics(const std::set<const Ribonucleotide*>& modifications,
                                                       const std::map<std::string, std::vector<double>>& explicit_mz,
                                                       bool negative_mode,
                                                       double tolerance,
                                                       bool tolerance_ppm,
                                                       double min_relative_intensity,
                                                       Int max_charge):
    tolerance_(tolerance),
    tolerance_ppm_(tolerance_ppm),
    min_relative_intensity_(min_relative_intensity)
{
  if (! std::isfinite(tolerance) || tolerance <= 0 || ! std::isfinite(min_relative_intensity) || min_relative_intensity < 0
      || min_relative_intensity > 1 || max_charge < 1 || max_charge > 6)
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid RNA diagnostic tolerance or intensity threshold");
  }
  const Int polarity = negative_mode ? -1 : 1;
  const double phosphate = EmpiricalFormula("HPO3").getMonoWeight();
  for (const auto* residue : modifications)
  {
    if (! residue || residue->getTermSpecificity() != Ribonucleotide::ANYWHERE) { continue; }
    Rule rule;
    rule.code = residue->getCode();
    const auto fixed = explicit_mz.find(rule.code);
    if (fixed != explicit_mz.end())
    {
      for (double mz : fixed->second)
      {
        rule.ions.push_back({mz, 1});
      }
      rule.labile = true;
    }
    else
    {
      // The terminal-state model also includes c+P/d+P single-nucleotide ions.
      // Use NASequence chemistry rather than a second table of residue masses.
      const NASequence single({residue}, nullptr, nullptr);
      for (Int charge = 1; charge <= max_charge; ++charge)
      {
        for (auto type : {NASequence::AIon, NASequence::BIon, NASequence::CIon, NASequence::DIon})
        {
          const double mz = single.getMonoWeight(type, polarity * charge) / charge;
          rule.ions.push_back({mz, charge});
          if (type == NASequence::CIon || type == NASequence::DIon) { rule.ions.push_back({mz + phosphate / charge, charge}); }
        }
      }
    }
    if (rule.ions.empty() || std::any_of(rule.ions.begin(), rule.ions.end(), [](const auto& ion) { return ! std::isfinite(ion.mz) || ion.mz <= 0; }))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RNA diagnostic m/z values must be finite and positive");
    }
    rules_.push_back(std::move(rule));
  }
}

RNAModificationDiagnostics::Evidence RNAModificationDiagnostics::evaluate(const MSSpectrum& raw) const
{
  Evidence evidence;
  for (const auto& rule : rules_)
  {
    evidence.tested.insert(rule.code);
  }
  double base_peak = 0;
  for (const auto& peak : raw)
  {
    if (std::isfinite(peak.getIntensity())) { base_peak = std::max(base_peak, double(peak.getIntensity())); }
  }
  if (base_peak <= 0) { return evidence; }
  std::vector<std::pair<double, double>> peaks;
  for (const auto& peak : raw)
  {
    if (std::isfinite(peak.getMZ()) && std::isfinite(peak.getIntensity()) && peak.getIntensity() > 0
        && peak.getIntensity() >= min_relative_intensity_ * base_peak)
    {
      peaks.emplace_back(peak.getMZ(), peak.getIntensity());
    }
  }
  std::sort(peaks.begin(), peaks.end());
  Int precursor_charge = 0;
  for (const auto& precursor : raw.getPrecursors())
  {
    precursor_charge = std::max(precursor_charge, std::abs(precursor.getCharge()));
  }
  for (const auto& rule : rules_)
  {
    for (const auto& ion : rule.ions)
    {
      if (precursor_charge && ion.charge > precursor_charge) { continue; }
      const double mz = ion.mz;
      const double delta = tolerance_ * (tolerance_ppm_ ? mz * 1e-6 : 1.0);
      auto found = std::lower_bound(peaks.begin(), peaks.end(), mz - delta, [](const auto& peak, double mass) { return peak.first < mass; });
      for (; found != peaks.end() && found->first <= mz + delta; ++found)
      {
        double isotope_match = 0;
        if (ion.charge > 1)
        {
          // A half-Dalton or smaller raw isotope doublet establishes charge;
          // a lone coincident peak must not relax modification permissions.
          const double partner_mz = found->first + Constants::C13C12_MASSDIFF_U / ion.charge;
          const double partner_delta = tolerance_ * (tolerance_ppm_ ? partner_mz * 1e-6 : 1.0);
          auto partner = std::lower_bound(peaks.begin(), peaks.end(), partner_mz - partner_delta,
                                          [](const auto& peak, double mass) { return peak.first < mass; });
          if (partner == found) { ++partner; }
          if (partner == peaks.end() || partner->first > partner_mz + partner_delta) { continue; }
          isotope_match = partner->first;
        }
        evidence.allowed.insert(rule.code);
        if (rule.labile) { evidence.labile.insert(rule.code); }
        evidence.matched_mz.emplace(rule.code, found->first);
        evidence.matched_charge.emplace(rule.code, ion.charge);
        if (isotope_match) { evidence.isotope_mz.emplace(rule.code, isotope_match); }
        break;
      }
      if (evidence.allowed.contains(rule.code)) { break; }
    }
  }
  return evidence;
}
} // namespace OpenMS

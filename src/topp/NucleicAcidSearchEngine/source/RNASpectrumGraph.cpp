// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "RNASpectrumGraph.h"

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

namespace OpenMS
{
RNASpectrumGraph::RNASpectrumGraph(const std::set<const Ribonucleotide*>& residues, const Options& options): options_(options)
{
  if (! std::isfinite(options.tolerance) || options.tolerance <= 0 || options.min_tag_length == 0 || options.max_tag_length < options.min_tag_length
      || options.max_tag_length > 32 || options.max_step < 1 || options.max_step > 2 || options.max_peaks < 2
      || options.max_peaks > std::numeric_limits<UInt>::max() || options.max_charge == 0 || options.max_edges == 0 || options.max_query_edges == 0
      || options.max_query_states == 0)
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid RNA graph tolerance, tag length or resource bound");
  }
  auto alphabet = std::make_shared<Alphabet>();
  if (options.tolerance_ppm && options.tolerance >= 1e6)
  {
    fallback_reason_ = "unsupported_tolerance";
    alphabet_ = std::move(alphabet);
    return;
  }
  std::map<std::string, double> masses;
  const double backbone = EmpiricalFormula("H-1PO2").getMonoWeight();
  for (const auto* residue : residues)
  {
    if (! residue || residue->getTermSpecificity() != Ribonucleotide::ANYWHERE || residue->getCode().empty() || residue->getCode().back() == '*')
    {
      continue; // Linkage-dependent increments require another graph model.
    }
    // Reuse the cached residue masses used by NucleicAcidSpectrumGenerator.
    // NASequence separately supplies formula masses for precursor validation.
    const double mass = residue->getMonoMass() + backbone;
    if (std::isfinite(mass) && mass > 0) { masses.emplace(residue->getCode(), mass); }
  }
  if (masses.empty()) { fallback_reason_ = "empty_alphabet"; }
  else if (masses.size() > (options.max_step == 2 ? 1023 : 65535)) { fallback_reason_ = "alphabet_limit"; }
  else
  {
    // Integer IDs are stable by code, and never collapse isobaric labels.
    for (const auto& [code, mass] : masses)
    {
      const UInt id = static_cast<UInt>(alphabet->code_ids.size());
      alphabet->code_ids.emplace(code, id);
      alphabet->transitions.push_back({mass, {id, 0}, 1});
    }
    if (options.max_step == 2)
    {
      const Size count = alphabet->transitions.size();
      alphabet->transitions.reserve(count + count * count);
      for (Size i = 0; i < count; ++i)
      {
        for (Size j = 0; j < count; ++j)
        {
          alphabet->transitions.push_back({alphabet->transitions[i].mass + alphabet->transitions[j].mass,
                                           {alphabet->transitions[i].residues[0], alphabet->transitions[j].residues[0]},
                                           2});
        }
      }
    }
    std::sort(alphabet->transitions.begin(), alphabet->transitions.end(),
              [](const auto& a, const auto& b) { return std::tie(a.mass, a.length, a.residues) < std::tie(b.mass, b.length, b.residues); });
  }
  alphabet_ = std::move(alphabet);
}

RNASpectrumGraph::Result
RNASpectrumGraph::build(const MSSpectrum& spectrum, const std::set<const Ribonucleotide*>& alphabet, const Options& options, Size max_charge)
{
  return RNASpectrumGraph(alphabet, options).build(spectrum, max_charge);
}

const MSSpectrum::IntegerDataArray* RNASpectrumGraph::legacyCharges(const MSSpectrum& spectrum)
{
  for (const auto& array : spectrum.getIntegerDataArrays())
  {
    if (array.getName() == "NASE:legacy_fragment_charge" && array.size() == spectrum.size()) { return &array; }
  }
  return nullptr;
}

std::vector<Size> RNASpectrumGraph::selectPeakIndices(const MSSpectrum& spectrum, Size max_peaks)
{
  std::vector<Size> peaks;
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    if (std::isfinite(spectrum[i].getMZ()) && spectrum[i].getMZ() > 0 && std::isfinite(spectrum[i].getIntensity()) && spectrum[i].getIntensity() > 0)
    {
      peaks.push_back(i);
    }
  }
  std::sort(peaks.begin(), peaks.end(), [&spectrum](Size a, Size b) {
    if (spectrum[a].getIntensity() != spectrum[b].getIntensity()) { return spectrum[a].getIntensity() > spectrum[b].getIntensity(); }
    return std::tuple(spectrum[a].getMZ(), a) < std::tuple(spectrum[b].getMZ(), b);
  });
  if (peaks.size() > max_peaks) { peaks.resize(max_peaks); }
  std::sort(peaks.begin(), peaks.end(),
            [&spectrum](Size a, Size b) { return std::tuple(spectrum[a].getMZ(), a) < std::tuple(spectrum[b].getMZ(), b); });
  return peaks;
}

RNASpectrumGraph::Result RNASpectrumGraph::build(const MSSpectrum& spectrum, Size max_charge) const
{
  if (max_charge == 0) { throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RNA graph fragment charge must be positive"); }
  Result result;
  result.alphabet = alphabet_;
  result.options = options_;
  result.fallback_reason = fallback_reason_;
  if (! result.fallback_reason.empty()) { return result; }
  result.peak_indices = selectPeakIndices(spectrum, options_.max_peaks);
  const auto& peaks = result.peak_indices;
  result.peaks = peaks.size();
  const auto& transitions = alphabet_->transitions;
  Size built_edges = 0;
  for (Size charge = 1; charge <= std::min(max_charge, options_.max_charge); ++charge)
  {
    Layer layer;
    layer.charge = charge;
    layer.outgoing.resize(peaks.size());
    Size layer_edges = 0;
    for (Size i = 0; i < peaks.size(); ++i)
    {
      for (Size j = i + 1; j < peaks.size(); ++j)
      {
        const double mz1 = spectrum[peaks[i]].getMZ(), mz2 = spectrum[peaks[j]].getMZ();
        const double mass = (mz2 - mz1) * charge;
        // Both endpoints contribute uncertainty; the charge carrier cancels.
        const double fraction = options_.tolerance * 1e-6;
        const double tolerance = charge * (options_.tolerance_ppm ? (mz1 + mz2) * fraction / (1.0 - fraction) : 2.0 * options_.tolerance);
        if (mass > transitions.back().mass + tolerance) { break; }
        const auto first = std::lower_bound(transitions.begin(), transitions.end(), mass - tolerance,
                                            [](const auto& transition, double value) { return transition.mass < value; });
        const auto past = std::upper_bound(first, transitions.end(), mass + tolerance,
                                           [](double value, const auto& transition) { return value < transition.mass; });
        for (UInt length = 1; length <= options_.max_step; ++length)
        {
          if (std::none_of(first, past, [length](const auto& transition) { return transition.length == length; })) { continue; }
          if (++built_edges > options_.max_edges)
          {
            result.layers.clear();
            result.peak_indices.clear();
            result.edges = built_edges;
            result.vertices = 0;
            result.fallback_reason = "edge_limit";
            return result;
          }
          ++layer_edges;
          layer.outgoing[i].push_back(
            {static_cast<UInt>(j), static_cast<UInt>(first - transitions.begin()), static_cast<UInt>(past - transitions.begin()), length});
        }
      }
    }
    bool layer_has_tag = false;
    {
      const Size width = options_.max_tag_length + 1;
      std::vector<unsigned char> reachable(peaks.size() * width, 0);
      for (Size node = 0; node < peaks.size() && ! layer_has_tag; ++node)
      {
        reachable[node * width] = 1;
        for (Size length = 0; length < width && ! layer_has_tag; ++length)
        {
          if (! reachable[node * width + length]) { continue; }
          for (const auto& edge : layer.outgoing[node])
          {
            const Size next = length + edge.length;
            if (next >= width) { continue; }
            reachable[edge.to * width + next] = 1;
            if (next >= options_.min_tag_length)
            {
              layer_has_tag = true;
              break;
            }
          }
        }
      }
    }
    if (layer_has_tag)
    {
      result.vertices += peaks.size();
      result.edges += layer_edges;
      result.layers.push_back(std::move(layer));
    }
  }
  if (result.layers.empty())
  {
    result.fallback_reason = "no_tags";
    result.layers.clear();
    result.peak_indices.clear();
  }
  return result;
}
} // namespace OpenMS

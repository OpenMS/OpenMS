// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "RNATagCandidateFilter.h"

#include <OpenMS/CONCEPT/Constants.h>
#include <algorithm>
#include <limits>
#include <charconv>
#include <cmath>
#include <map>
#include <string_view>
#include <system_error>
#include <type_traits>

namespace OpenMS
{
namespace
{
  bool prefixIon(std::string_view ion)
  {
    return ion == "a" || ion == "b" || ion == "c" || ion == "d";
  }
  bool suffixIon(std::string_view ion)
  {
    return ion == "w" || ion == "x" || ion == "y" || ion == "z";
  }

  std::string fallback(const NASequence& sequence, const RNASpectrumGraph::Result& graph, const std::set<std::string>& ions)
  {
    if (! graph.fallback_reason.empty()) { return graph.fallback_reason; }
    if (! graph.alphabet || graph.layers.empty() || graph.peak_indices.empty()) { return "no_tags"; }
    if (sequence.size() <= graph.options.min_tag_length + 1) { return "short_candidate"; }
    if (std::none_of(ions.begin(), ions.end(), [](const auto& ion) { return prefixIon(ion) || suffixIon(ion); }))
    {
      return "no_supported_ion_series";
    }
    for (const auto& residue : sequence)
    {
      if (! graph.alphabet->code_ids.contains(residue.getCode())) { return "unsupported_residue"; }
    }
    return {};
  }

  std::vector<UInt> residueIDs(const NASequence& sequence, const RNASpectrumGraph::Alphabet& alphabet)
  {
    std::vector<UInt> result;
    result.reserve(sequence.size());
    for (const auto& residue : sequence)
    {
      result.push_back(alphabet.code_ids.at(residue.getCode()));
    }
    return result;
  }

  /// One state bound for the prefilter and both matchers: a candidate that a
  /// matcher keeps through "query_state_limit" must not be rejected earlier.
  bool exceedsQueryStates(const RNASpectrumGraph::Result& graph, Size sequence_length)
  {
    const Size nodes = graph.peak_indices.size(), width = graph.options.max_tag_length + 1;
    if (nodes == 0) { return false; }
    return nodes > graph.options.max_query_states / width || sequence_length > graph.options.max_query_states / (nodes * width);
  }

  enum class StepMatch
  {
    MATCH,
    NONE,
    LIMIT
  };

  bool containsResidue(UInt code, UInt residue)
  {
    return code == residue;
  }
  bool containsResidue(const std::vector<UInt>& choices, UInt residue)
  {
    return std::binary_search(choices.begin(), choices.end(), residue);
  }

  template<typename Codes>
  StepMatch matchStep(const RNASpectrumGraph::Alphabet& alphabet,
                      const RNASpectrumGraph::Edge& edge,
                      const Codes& codes,
                      Size position,
                      bool suffix,
                      Size& work,
                      Size limit)
  {
    bool matched = false;
    for (Size i = edge.first_transition; i < edge.past_transition; ++i)
    {
      if (++work > limit) { return StepMatch::LIMIT; }
      const auto& transition = alphabet.transitions[i];
      if (transition.length != edge.length) { continue; }
      bool same = true;
      for (Size offset = 0; offset < edge.length; ++offset)
      {
        const Size index = suffix ? codes.size() - position - offset - 1 : position + offset;
        if (! containsResidue(codes[index], transition.residues[offset]))
        {
          same = false;
          break;
        }
      }
      if (same)
      {
        if constexpr (std::is_same_v<typename Codes::value_type, UInt>) { return StepMatch::MATCH; }
        matched = true;
      }
    }
    // Domain queries must bound every singleton completion's work, including
    // isobars found later in the transition range. Otherwise a partial query
    // could reject a completion that the reference keeps through query fallback.
    return matched ? StepMatch::MATCH : StepMatch::NONE;
  }
  template<typename Codes>
  bool queryCompatible(const Codes& codes, const RNASpectrumGraph::Result& graph, const std::set<std::string>& ions, Size max_charge)
  {
    if (max_charge && std::none_of(graph.layers.begin(), graph.layers.end(), [max_charge](const auto& layer) { return layer.charge <= max_charge; }))
    {
      return true;
    }
    const auto& options = graph.options;
    const Size width = options.max_tag_length + 1;
    const Size nodes = graph.peak_indices.size();
    if (exceedsQueryStates(graph, codes.size())) { return true; } // The matcher would fall back.
    Size work = 0;
    for (bool suffix : {false, true})
    {
      if (std::none_of(ions.begin(), ions.end(), [suffix](const auto& ion) { return suffix ? suffixIon(ion) : prefixIon(ion); })) { continue; }
      for (const auto& layer : graph.layers)
      {
        if (max_charge && layer.charge > max_charge) { continue; }
        for (Size start = 1; start + options.min_tag_length < codes.size(); ++start)
        {
          // Most observed steps (plus one; zero is unreachable) per node and residue span.
          std::vector<unsigned char> steps(nodes * width, 0);
          for (Size node = 0; node < nodes; ++node)
          {
            steps[node * width] = 1; // Any observed ion may start an internal tag.
            for (Size length = 0; length < width; ++length)
            {
              const Size taken = steps[node * width + length];
              if (! taken) { continue; }
              for (const auto& edge : layer.outgoing[node])
              {
                const Size next = length + edge.length;
                if (next >= width || start + next >= codes.size()) { continue; }
                if (++work > options.max_query_edges)
                {
                  return true; // Incomplete query fails open.
                }
                const auto matched = matchStep(*graph.alphabet, edge, codes, start + length, suffix, work, options.max_query_edges);
                if (matched == StepMatch::LIMIT) { return true; }
                if (matched == StepMatch::NONE) { continue; }
                if (taken >= options.min_tag_length) { return true; }
                auto& best = steps[edge.to * width + next];
                best = std::max<unsigned char>(best, taken + 1);
              }
            }
          }
        }
      }
    }
    return false;
  }
} // namespace

bool RNATagCandidateFilter::isCompatible(const NASequence& candidate,
                                         const RNASpectrumGraph::Result& graph,
                                         const std::set<std::string>& ions,
                                         Size max_charge)
{
  if (! fallback(candidate, graph, ions).empty()) { return true; }
  return queryCompatible(residueIDs(candidate, *graph.alphabet), graph, ions, max_charge);
}

bool RNATagCandidateFilter::isCompatible(const std::vector<std::vector<const Ribonucleotide*>>& alternatives,
                                         const RNASpectrumGraph::Result& graph,
                                         const std::set<std::string>& ions,
                                         Size max_charge)
{
  if (! graph.fallback_reason.empty() || ! graph.alphabet || graph.layers.empty() || graph.peak_indices.empty()
      || alternatives.size() <= graph.options.min_tag_length + 1
      || std::none_of(ions.begin(), ions.end(), [](const auto& ion) { return prefixIon(ion) || suffixIon(ion); }))
  {
    return true;
  }
  std::vector<std::vector<UInt>> codes;
  codes.reserve(alternatives.size());
  for (const auto& choices : alternatives)
  {
    if (choices.empty()) { return false; }
    auto& ids = codes.emplace_back();
    for (const auto* residue : choices)
    {
      if (! residue) { return true; }
      const auto found = graph.alphabet->code_ids.find(residue->getCode());
      if (found == graph.alphabet->code_ids.end()) { return true; }
      ids.push_back(found->second);
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  }
  return queryCompatible(codes, graph, ions, max_charge);
}


RNATagCandidateFilter::PreparedQuery RNATagCandidateFilter::prepare(const NASequence& candidate,
                                                                    std::shared_ptr<const RNASpectrumGraph::Alphabet> alphabet,
                                                                    const MSSpectrum& theoretical,
                                                                    const std::set<std::string>& ions)
{
  PreparedQuery query;
  query.alphabet = std::move(alphabet);
  query.sequence_length = candidate.size();
  query.supported_series = std::any_of(ions.begin(), ions.end(), [](const auto& ion) { return prefixIon(ion) || suffixIon(ion); });
  query.supported_residues = bool(query.alphabet);
  if (query.alphabet)
  {
    query.codes.reserve(candidate.size());
    for (const auto& residue : candidate)
    {
      const auto found = query.alphabet->code_ids.find(residue.getCode());
      if (found == query.alphabet->code_ids.end())
      {
        query.supported_residues = false;
        break;
      }
      query.codes.push_back(found->second);
    }
  }
  if (! query.supported_series || ! query.supported_residues) { return query; }
  const MSSpectrum::StringDataArray* names = nullptr;
  const MSSpectrum::IntegerDataArray* charges = nullptr;
  for (const auto& array : theoretical.getStringDataArrays())
  {
    if (array.getName() == Constants::UserParam::IonNames) { names = &array; }
  }
  for (const auto& array : theoretical.getIntegerDataArrays())
  {
    if (array.getName() == "Charges") { charges = &array; }
  }
  query.valid_metadata = names && charges && names->size() == theoretical.size() && charges->size() == theoretical.size();
  if (! query.valid_metadata) { return query; }
  for (Int charge : *charges)
  {
    query.theoretical_charges.push_back(static_cast<Size>(std::abs(charge)));
  }
  std::sort(query.theoretical_charges.begin(), query.theoretical_charges.end());
  query.theoretical_charges.erase(std::unique(query.theoretical_charges.begin(), query.theoretical_charges.end()), query.theoretical_charges.end());
  std::map<std::pair<char, Size>, std::vector<PreparedQuery::Ion>> grouped;
  for (Size i = 0; i < theoretical.size(); ++i)
  {
    const auto& name = (*names)[i];
    if (name.size() < 2) { continue; }
    const std::string series(1, name[0]);
    if ((! prefixIon(series) && ! suffixIon(series)) || ! ions.contains(series)) { continue; }
    Size position = 0;
    const auto parsed = std::from_chars(name.data() + 1, name.data() + name.size(), position);
    if (parsed.ec != std::errc {} || parsed.ptr != name.data() + name.size() || position == 0 || position >= candidate.size()) { continue; }
    grouped[{name[0], static_cast<Size>(std::abs((*charges)[i]))}].push_back({position, theoretical[i].getMZ()});
  }
  for (auto& [key, peaks] : grouped)
  {
    const bool sorted_finite = std::all_of(peaks.begin(), peaks.end(), [](const auto& peak) { return std::isfinite(peak.mz); })
                               && std::is_sorted(peaks.begin(), peaks.end(), [](const auto& a, const auto& b) { return a.mz < b.mz; });
    query.groups.push_back({key.first, key.second, std::move(peaks), sorted_finite});
  }
  return query;
}

RNATagCandidateFilter::Evidence RNATagCandidateFilter::match(const NASequence& candidate,
                                                             const RNASpectrumGraph::Result& graph,
                                                             const MSSpectrum& experimental,
                                                             const MSSpectrum& theoretical,
                                                             const std::set<std::string>& ions,
                                                             double tolerance,
                                                             bool tolerance_ppm)
{
  return match(prepare(candidate, graph.alphabet, theoretical, ions), graph, experimental, tolerance, tolerance_ppm);
}

RNATagCandidateFilter::Evidence RNATagCandidateFilter::match(const PreparedQuery& query,
                                                             const RNASpectrumGraph::Result& graph,
                                                             const MSSpectrum& experimental,
                                                             double tolerance,
                                                             bool tolerance_ppm)
{
  Evidence evidence;
  const auto fallbackTo = [&evidence](const std::string& reason) {
    evidence.fallback_reason = reason;
    return evidence;
  };
  if (! graph.fallback_reason.empty()) { return fallbackTo(graph.fallback_reason); }
  if (! graph.alphabet || graph.layers.empty() || graph.peak_indices.empty()) { return fallbackTo("no_tags"); }
  if (query.sequence_length <= graph.options.min_tag_length + 1) { return fallbackTo("short_candidate"); }
  if (! query.supported_series) { return fallbackTo("no_supported_ion_series"); }
  if (! query.supported_residues) { return fallbackTo("unsupported_residue"); }
  const Size width = graph.options.max_tag_length + 1;
  if (exceedsQueryStates(graph, query.sequence_length)) { return fallbackTo("query_state_limit"); }
  if (std::any_of(graph.peak_indices.begin(), graph.peak_indices.end(), [&experimental](Size index) { return index >= experimental.size(); }))
  {
    return fallbackTo("invalid_peak_index");
  }
  if (! query.valid_metadata) { return fallbackTo("missing_theoretical_annotations"); }
  if (std::none_of(graph.layers.begin(), graph.layers.end(), [&query](const auto& layer) {
        return std::binary_search(query.theoretical_charges.begin(), query.theoretical_charges.end(), layer.charge);
      }))
  {
    return fallbackTo("no_tags_for_charge");
  }
  if (query.alphabet != graph.alphabet) { return fallbackTo("alphabet_mismatch"); }

  // A group only needs states for observed ion matches. Sparse state preserves
  // the reference's traversal order and query-work accounting exactly.
  std::vector<std::pair<Size, UInt>> matches;
  std::vector<Size> row_offsets;
  // Most observed steps (plus one; zero is unreachable) per match and residue span.
  std::vector<unsigned char> steps;
  Size work = 0;
  for (const auto& group : query.groups)
  {
    const auto layer_it
      = std::find_if(graph.layers.begin(), graph.layers.end(), [&group](const auto& layer) { return layer.charge == group.charge; });
    if (layer_it == graph.layers.end()) { continue; }
    const auto& layer = *layer_it;
    auto first_ion = group.ions.begin(), past_ion = group.ions.end();
    if (group.sorted_finite_mz && std::isfinite(tolerance) && tolerance >= 0)
    {
      const double minimum = experimental[graph.peak_indices.front()].getMZ();
      const double maximum = experimental[graph.peak_indices.back()].getMZ();
      double low = minimum - tolerance, high = maximum + tolerance;
      const double fraction = tolerance * 1e-6;
      if (! tolerance_ppm || fraction < 1.0)
      {
        if (tolerance_ppm)
        {
          low = minimum / (1.0 + fraction);
          high = maximum / (1.0 - fraction);
        }
        // Outward slack keeps endpoint rounding from dropping an exact match.
        const double slack = 32 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(low), std::abs(high)});
        first_ion = std::lower_bound(first_ion, past_ion, low - slack, [](const auto& ion, double value) { return ion.mz < value; });
        past_ion = std::upper_bound(first_ion, past_ion, high + slack, [](double value, const auto& ion) { return value < ion.mz; });
      }
    }
    matches.clear();
    for (auto ion = first_ion; ion != past_ion; ++ion)
    {
      const double delta = tolerance * (tolerance_ppm ? ion->mz * 1e-6 : 1.0);
      auto first = std::lower_bound(graph.peak_indices.begin(), graph.peak_indices.end(), ion->mz - delta,
                                    [&experimental](Size index, double value) { return experimental[index].getMZ() < value; });
      for (auto it = first; it != graph.peak_indices.end() && experimental[*it].getMZ() <= ion->mz + delta; ++it)
      {
        matches.emplace_back(ion->position, static_cast<UInt>(it - graph.peak_indices.begin()));
      }
    }
    if (matches.empty()) { continue; }
    std::sort(matches.begin(), matches.end());
    matches.erase(std::unique(matches.begin(), matches.end()), matches.end());
    row_offsets.assign(query.sequence_length + 1, 0);
    for (const auto& [position, node] : matches)
    {
      ++row_offsets[position + 1];
    }
    for (Size i = 1; i < row_offsets.size(); ++i)
    {
      row_offsets[i] += row_offsets[i - 1];
    }
    steps.assign(matches.size() * width, 0);
    const bool suffix = group.series >= 'w';
    for (Size position = 1; position < query.sequence_length; ++position)
    {
      for (Size index = row_offsets[position]; index < row_offsets[position + 1]; ++index)
      {
        const UInt node = matches[index].second;
        steps[index * width] = 1;
        for (const auto& edge : layer.outgoing[node])
        {
          if (++work > graph.options.max_query_edges)
          {
            if (! evidence.supported) { evidence.fallback_reason = "query_limit"; }
            return evidence;
          }
          const Size next_position = position + edge.length;
          if (next_position >= query.sequence_length) { continue; }
          const auto begin = matches.begin() + row_offsets[next_position], end = matches.begin() + row_offsets[next_position + 1];
          const auto next = std::lower_bound(begin, end, edge.to, [](const auto& match, UInt value) { return match.second < value; });
          if (next == end || next->second != edge.to) { continue; }
          const auto matched = matchStep(*graph.alphabet, edge, query.codes, position, suffix, work, graph.options.max_query_edges);
          if (matched == StepMatch::LIMIT)
          {
            if (! evidence.supported) { evidence.fallback_reason = "query_limit"; }
            return evidence;
          }
          if (matched == StepMatch::NONE) { continue; }
          const Size next_index = next - matches.begin();
          for (Size length = 0; length + edge.length < width; ++length)
          {
            const Size taken = steps[index * width + length];
            if (! taken) { continue; }
            const Size next_length = length + edge.length;
            auto& best = steps[next_index * width + next_length];
            best = std::max<unsigned char>(best, taken + 1);
            if (taken >= graph.options.min_tag_length && next_length > evidence.length)
            {
              evidence.supported = true;
              evidence.length = next_length;
              evidence.steps = taken;
              evidence.charge = layer.charge;
              evidence.ion_series.assign(1, group.series);
              evidence.start_ion = next_position - next_length;
            }
          }
          if (evidence.length == graph.options.max_tag_length) { return evidence; }
        }
      }
    }
  }
  return evidence;
}

RNATagCandidateFilter::Evidence RNATagCandidateFilter::matchReference(const NASequence& candidate,
                                                                      const RNASpectrumGraph::Result& graph,
                                                                      const MSSpectrum& experimental,
                                                                      const MSSpectrum& theoretical,
                                                                      const std::set<std::string>& ions,
                                                                      double tolerance,
                                                                      bool tolerance_ppm)
{
  Evidence evidence;
  evidence.fallback_reason = fallback(candidate, graph, ions);
  if (! evidence.fallback_reason.empty()) { return evidence; }
  const Size nodes = graph.peak_indices.size(), width = graph.options.max_tag_length + 1;
  if (exceedsQueryStates(graph, candidate.size()))
  {
    evidence.fallback_reason = "query_state_limit";
    return evidence;
  }
  if (std::any_of(graph.peak_indices.begin(), graph.peak_indices.end(), [&experimental](Size index) { return index >= experimental.size(); }))
  {
    evidence.fallback_reason = "invalid_peak_index";
    return evidence;
  }
  const MSSpectrum::StringDataArray* names = nullptr;
  const MSSpectrum::IntegerDataArray* charges = nullptr;
  for (const auto& array : theoretical.getStringDataArrays())
  {
    if (array.getName() == Constants::UserParam::IonNames) { names = &array; }
  }
  for (const auto& array : theoretical.getIntegerDataArrays())
  {
    if (array.getName() == "Charges") { charges = &array; }
  }
  if (! names || ! charges || names->size() != theoretical.size() || charges->size() != theoretical.size())
  {
    evidence.fallback_reason = "missing_theoretical_annotations";
    return evidence;
  }
  std::set<Size> theoretical_charges;
  for (Int charge : *charges)
  {
    theoretical_charges.insert(static_cast<Size>(std::abs(charge)));
  }
  if (std::none_of(graph.layers.begin(), graph.layers.end(),
                   [&theoretical_charges](const auto& layer) { return theoretical_charges.contains(layer.charge); }))
  {
    evidence.fallback_reason = "no_tags_for_charge";
    return evidence;
  }

  // Find observed nodes for the generator's actual ion masses. This reuses
  // terminal offsets and carrier masses; it never joins different ion series.
  using IonKey = std::pair<std::string, Size>;
  std::map<IonKey, std::vector<std::vector<UInt>>> matches;
  for (Size i = 0; i < theoretical.size(); ++i)
  {
    const std::string& name = (*names)[i];
    if (name.size() < 2) { continue; }
    const std::string series(1, name[0]);
    if ((! prefixIon(series) && ! suffixIon(series)) || ! ions.contains(series)) { continue; }
    Size position = 0;
    const auto parsed = std::from_chars(name.data() + 1, name.data() + name.size(), position);
    if (parsed.ec != std::errc {} || parsed.ptr != name.data() + name.size() || position == 0 || position >= candidate.size()) { continue; }
    const Size charge = std::abs((*charges)[i]);
    if (std::none_of(graph.layers.begin(), graph.layers.end(), [charge](const auto& layer) { return layer.charge == charge; })) { continue; }
    const double mass = theoretical[i].getMZ();
    const double delta = tolerance * (tolerance_ppm ? mass * 1e-6 : 1.0);
    auto first = std::lower_bound(graph.peak_indices.begin(), graph.peak_indices.end(), mass - delta,
                                  [&experimental](Size index, double value) { return experimental[index].getMZ() < value; });
    for (auto it = first; it != graph.peak_indices.end() && experimental[*it].getMZ() <= mass + delta; ++it)
    {
      auto& found = matches[{series, charge}];
      if (found.empty()) { found.resize(candidate.size()); }
      found[position].push_back(static_cast<UInt>(it - graph.peak_indices.begin()));
    }
  }
  const auto codes = residueIDs(candidate, *graph.alphabet);
  Size work = 0;
  const Size stride = nodes * width;
  for (auto& [key, found] : matches)
  {
    const auto charge = key.second;
    const auto& layer = *std::find_if(graph.layers.begin(), graph.layers.end(), [charge](const auto& item) { return item.charge == charge; });
    for (auto& row : found)
    {
      std::sort(row.begin(), row.end());
      row.erase(std::unique(row.begin(), row.end()), row.end());
    }
    // Most observed steps (plus one; zero is unreachable) per position, node and residue span.
    std::vector<unsigned char> steps(candidate.size() * stride, 0);
    for (Size position = 1; position < candidate.size(); ++position)
    {
      for (UInt node : found[position])
      {
        steps[position * stride + node * width] = 1;
        for (const auto& edge : layer.outgoing[node])
        {
          if (++work > graph.options.max_query_edges)
          {
            if (! evidence.supported) { evidence.fallback_reason = "query_limit"; }
            return evidence;
          }
          const Size next_position = position + edge.length;
          if (next_position >= candidate.size() || ! std::binary_search(found[next_position].begin(), found[next_position].end(), edge.to))
          {
            continue;
          }
          const auto matched = matchStep(*graph.alphabet, edge, codes, position, suffixIon(key.first), work, graph.options.max_query_edges);
          if (matched == StepMatch::LIMIT)
          {
            if (! evidence.supported) { evidence.fallback_reason = "query_limit"; }
            return evidence;
          }
          if (matched == StepMatch::NONE) { continue; }
          for (Size length = 0; length + edge.length < width; ++length)
          {
            const Size taken = steps[position * stride + node * width + length];
            if (! taken) { continue; }
            const Size next_length = length + edge.length;
            auto& best = steps[next_position * stride + edge.to * width + next_length];
            best = std::max<unsigned char>(best, taken + 1);
            if (taken >= graph.options.min_tag_length && next_length > evidence.length)
            {
              evidence.supported = true;
              evidence.length = next_length;
              evidence.steps = taken;
              evidence.charge = layer.charge;
              evidence.ion_series = key.first;
              evidence.start_ion = next_position - next_length;
            }
          }
        }
      }
    }
  }
  return evidence;
}
} // namespace OpenMS

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "RNARescoringFeatures.h"

#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>

namespace OpenMS
{
RNARescoringFeatures::RNARescoringFeatures(Options options) : options_(std::move(options))
{
  if (options_.min_charge == 0 || options_.min_charge > options_.max_charge)
  {
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Rescoring features need a charge range 0 < min_charge <= max_charge");
  }
  if (options_.route_score && ! options_.route)
  {
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The route_score column needs the route features");
  }
  for (const auto& series : options_.ion_series)
  {
    if (! (series == "a-B" || (series.size() == 1 && std::string("abcdwxyz").find(series) != std::string::npos)))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown ion series '" + series + "'");
    }
  }
  names_ = {"native_score", "delta_next", "delta_best", "abs_ppm_error", "isotope_offset"};
  for (Size charge = options_.min_charge; charge <= options_.max_charge; ++charge)
  {
    names_.push_back("charge_" + std::to_string(charge));
  }
  names_.insert(names_.end(), {"length", "modifications", "adduct"});
  for (const auto& series : options_.ion_series) { names_.push_back("frac_" + series); }
  names_.insert(names_.end(), {"matched_ions", "matched_intensity", "longest_series", "coverage"});
  if (options_.route)
  {
    if (options_.route_score) { names_.push_back("route_score"); }
    names_.insert(names_.end(), {"route_nodes", "route_intensity_ratio", "route_prefix_length", "route_suffix_length"});
  }
  if (options_.graph) { names_.push_back("graph_tag_length"); }
  // the unsupported count would be modifications - diagnostic_supported_mods
  if (options_.diagnostics) { names_.push_back("diagnostic_supported_mods"); }
  if (options_.flash) { names_.insert(names_.end(), {"flash_corrected", "flash_cosine", "flash_snr"}); }
}

std::optional<std::pair<std::string, Size>> RNARescoringFeatures::parseAnnotation(const std::string& annotation)
{
  if (annotation.empty() || std::string("abcdwxyz").find(annotation[0]) == std::string::npos) { return std::nullopt; }
  Size end = 1;
  while (end < annotation.size() && std::isdigit(static_cast<unsigned char>(annotation[end]))) { ++end; }
  if (end == 1) { return std::nullopt; }
  const std::string rest = annotation.substr(end);
  std::string series(1, annotation[0]);
  if (rest == "-B")
  {
    if (series != "a") { return std::nullopt; }
    series = "a-B";
  }
  else if (! rest.empty() && rest.find_first_not_of("+-") != std::string::npos)
  {
    return std::nullopt; // e.g. neutral losses other than a-B: not a backbone ladder rung
  }
  return std::make_pair(series, static_cast<Size>(std::stoul(annotation.substr(1, end - 1))));
}

bool RNARescoringFeatures::isPrefixSeries(const std::string& series)
{
  return series == "a-B" || series == "a" || series == "b" || series == "c" || series == "d";
}

void RNARescoringFeatures::appendIonFeatures_(const HitInput& hit, const SpectrumInput& spectrum,
                                              std::vector<double>& row) const
{
  const Size length = hit.sequence.size();
  const Size rungs = length > 1 ? length - 1 : 0;
  std::map<std::string, std::set<Size>> ordinals; // per configured series, in [1, length - 1]
  for (const auto& series : options_.ion_series) { ordinals[series]; }
  std::set<double> matched_peaks; // a peak matched by several ions counts once
  double matched_intensity = 0.0;
  if (hit.annotations)
  {
    for (const auto& annotation : *hit.annotations)
    {
      const auto parsed = parseAnnotation(annotation.annotation);
      if (! parsed || parsed->second == 0 || parsed->second > rungs) { continue; }
      auto found = ordinals.find(parsed->first);
      if (found == ordinals.end()) { continue; }
      found->second.insert(parsed->second);
      if (matched_peaks.insert(annotation.mz).second) { matched_intensity += annotation.intensity; }
    }
  }

  const double denominator = rungs > 0 ? double(rungs) : 1.0;
  Size matched_ions = 0, longest = 0;
  std::set<Size> covered; // bond after residue p, p in [1, length - 1]
  for (const auto& series : options_.ion_series)
  {
    const auto& matched = ordinals[series];
    row.push_back(matched.size() / denominator);
    matched_ions += matched.size();
    longest = std::max(longest, SpectrumAnnotator::computeLongestIonSeries(matched.begin(), matched.end(), 0,
                                                                            [](Size ordinal) { return ordinal; }));
    const bool prefix = isPrefixSeries(series);
    for (Size ordinal : matched) { covered.insert(prefix ? ordinal : length - ordinal); }
  }
  const double possible_ions = denominator * std::max<Size>(options_.ion_series.size(), 1);
  row.push_back(matched_ions / possible_ions);
  row.push_back(spectrum.total_intensity > 0.0 ? matched_intensity / spectrum.total_intensity : 0.0);
  row.push_back(longest / denominator);
  row.push_back(covered.size() / denominator);
}

std::vector<double> RNARescoringFeatures::computeRow(const std::vector<HitInput>& hits, Size index,
                                                    const SpectrumInput& spectrum) const
{
  OPENMS_PRECONDITION(index < hits.size(), "The hit must be one of the hits of the spectrum");
  double best = 0.0;
  for (const auto& hit : hits) { best = std::max(best, hit.native_score); }
  // The margin compares one hit with all others: linear in the number of hits.
  const std::string sequence = hits[index].sequence.toString();
  double next = 0.0;
  for (Size other = 0; other < hits.size(); ++other)
  {
    const double score = hits[other].native_score;
    if ((score <= hits[index].native_score) && (score > next) && (hits[other].sequence.toString() != sequence))
    {
      next = score;
    }
  }
  return buildRow_(hits[index], next, best, spectrum);
}

std::vector<std::vector<double>> RNARescoringFeatures::computeSpectrum(const std::vector<HitInput>& hits,
                                                                       const SpectrumInput& spectrum) const
{
  double best = 0.0;
  for (const auto& hit : hits) { best = std::max(best, hit.native_score); }
  std::vector<std::string> sequences;
  sequences.reserve(hits.size());
  for (const auto& hit : hits) { sequences.push_back(hit.sequence.toString()); }

  std::vector<std::vector<double>> rows;
  rows.reserve(hits.size());
  for (Size index = 0; index < hits.size(); ++index)
  {
    double next = 0.0;
    for (Size other = 0; other < hits.size(); ++other)
    {
      const double score = hits[other].native_score;
      if (sequences[other] != sequences[index] && score <= hits[index].native_score) { next = std::max(next, score); }
    }
    rows.push_back(buildRow_(hits[index], next, best, spectrum));
  }
  return rows;
}

std::vector<double> RNARescoringFeatures::buildRow_(const HitInput& hit, double next, double best,
                                                    const SpectrumInput& spectrum) const
{
  std::vector<double> row;
  row.reserve(names_.size());
  row.insert(row.end(), {hit.native_score, hit.native_score - next, hit.native_score - best,
                         std::fabs(hit.precursor_error_ppm), double(hit.isotope_offset)});
  const Size charge = static_cast<Size>(std::abs(hit.charge));
  for (Size z = options_.min_charge; z <= options_.max_charge; ++z) { row.push_back(charge == z ? 1.0 : 0.0); }

  Size modifications = 0, supported = 0;
  for (const Ribonucleotide& residue : hit.sequence)
  {
    if (! residue.isModified()) { continue; }
    ++modifications;
    if (spectrum.diagnostic_allowed.count(residue.getCode())) { ++supported; }
  }
  row.insert(row.end(), {double(hit.sequence.size()), double(modifications), hit.adduct ? 1.0 : 0.0});
  appendIonFeatures_(hit, spectrum, row);

  if (options_.route)
  {
    if (options_.route_score) { row.push_back(hit.route_score); }
    row.insert(row.end(), {double(hit.route_nodes),
                           hit.route_local_intensity > 0.0 ? hit.route_matched_intensity / hit.route_local_intensity : 0.0,
                           double(hit.route_prefix_length), double(hit.route_suffix_length)});
  }
  if (options_.graph) { row.push_back(double(hit.graph_tag_length)); }
  if (options_.diagnostics) { row.push_back(double(supported)); }
  if (options_.flash)
  {
    const bool corrected = spectrum.flash_cosine.has_value();
    row.insert(row.end(), {corrected ? 1.0 : 0.0, spectrum.flash_cosine.value_or(0.0), spectrum.flash_snr.value_or(0.0)});
  }
  OPENMS_POSTCONDITION(row.size() == names_.size(), "Every feature row must match the feature names");
  return row;
}
} // namespace OpenMS

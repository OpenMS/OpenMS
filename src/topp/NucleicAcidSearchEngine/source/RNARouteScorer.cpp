// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "RNARouteScorer.h"

#include <OpenMS/CHEMISTRY/SpectrumAnnotator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <tuple>

namespace OpenMS
{
RNARouteScorer::PreparedSpectrum RNARouteScorer::prepareSpectrum(const MSSpectrum& experimental, bool negative_mode, const Options& options)
{
  if (! std::isfinite(options.local_radius) || options.local_radius <= 0 || options.local_peaks == 0 || options.max_unknown_charge == 0)
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid RNA route background options");
  }
  PreparedSpectrum prepared;
  prepared.backgrounds.resize(experimental.size());
  const auto* charges = RNASpectrumGraph::legacyCharges(experimental);
  const double carrier = negative_mode ? -Constants::PROTON_MASS_U : Constants::PROTON_MASS_U;
  struct Node
  {
    double mass;
    Size peak, charge;
  };
  std::vector<Node> nodes;
  for (Size i = 0; i < experimental.size(); ++i)
  {
    if (! std::isfinite(experimental[i].getMZ()) || ! std::isfinite(experimental[i].getIntensity()) || experimental[i].getIntensity() <= 0)
    {
      continue;
    }
    const Size assigned = charges ? Size(std::abs(Int64((*charges)[i]))) : 0;
    const Size first = assigned ? assigned : 1, last = assigned ? assigned : options.max_unknown_charge;
    for (Size charge = first; charge <= last; ++charge)
    {
      const double mass = (experimental[i].getMZ() - carrier) * charge;
      if (mass > 0) { nodes.push_back({mass, i, charge}); }
    }
  }
  std::sort(nodes.begin(), nodes.end(),
            [](const auto& a, const auto& b) { return std::tie(a.mass, a.peak, a.charge) < std::tie(b.mass, b.peak, b.charge); });
  for (const auto& node : nodes)
  {
    std::vector<Size> local;
    const auto first = std::lower_bound(nodes.begin(), nodes.end(), node.mass - options.local_radius,
                                        [](const auto& value, double mass) { return value.mass < mass; });
    for (auto it = first; it != nodes.end() && it->mass <= node.mass + options.local_radius; ++it)
    {
      local.push_back(it->peak);
    }
    std::sort(local.begin(), local.end());
    local.erase(std::unique(local.begin(), local.end()), local.end());
    std::sort(local.begin(), local.end(),
              [&](Size a, Size b) { return std::tuple(-experimental[a].getIntensity(), a) < std::tuple(-experimental[b].getIntensity(), b); });
    if (local.size() > options.local_peaks) { local.resize(options.local_peaks); }
    if (std::find(local.begin(), local.end(), node.peak) == local.end()) { local.push_back(node.peak); }
    prepared.backgrounds[node.peak].emplace(node.charge, std::move(local));
  }
  return prepared;
}

RNARouteScorer::Result
RNARouteScorer::score(const RNATagCandidateFilter::PreparedQuery& query, const MSSpectrum& experimental, bool negative_mode, const Options& options)
{
  return score(query, experimental, prepareSpectrum(experimental, negative_mode, options), negative_mode, options);
}

RNARouteScorer::Result RNARouteScorer::score(const RNATagCandidateFilter::PreparedQuery& query,
                                             const MSSpectrum& experimental,
                                             const PreparedSpectrum& prepared,
                                             bool negative_mode,
                                             const Options& options)
{
  const auto routes = match(query, experimental, prepared, negative_mode, options);
  return score(std::vector<const PreparedRoutes*> {&routes}, experimental, prepared, negative_mode);
}

RNARouteScorer::PreparedRoutes RNARouteScorer::match(const RNATagCandidateFilter::PreparedQuery& query,
                                                     const MSSpectrum& experimental,
                                                     const PreparedSpectrum& prepared,
                                                     bool negative_mode,
                                                     const Options& options,
                                                     const std::string& scaffold)
{
  if (! std::isfinite(options.tolerance) || options.tolerance <= 0 || ! std::isfinite(options.local_radius) || options.local_radius <= 0
      || options.local_peaks == 0 || options.max_step == 0 || options.min_nodes == 0 || options.max_unknown_charge == 0
      || ! std::isfinite(options.mass_accuracy_tolerance) || options.mass_accuracy_tolerance < 0)
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid RNA route score options");
  }
  PreparedRoutes result;
  if (! query.valid_metadata || ! experimental.isSorted() || prepared.backgrounds.size() != experimental.size()) { return result; }
  std::map<char, std::vector<MatchedNode>> by_series;
  const double carrier = negative_mode ? -Constants::PROTON_MASS_U : Constants::PROTON_MASS_U;
  for (const auto& group : query.groups)
  {
    for (const auto& ion : group.ions)
    {
      const double delta = options.tolerance * (options.tolerance_ppm ? ion.mz * 1e-6 : 1.0);
      const auto first = experimental.MZBegin(ion.mz - delta);
      for (auto peak = first; peak != experimental.end() && peak->getMZ() <= ion.mz + delta; ++peak)
      {
        if (! std::isfinite(peak->getIntensity()) || peak->getIntensity() <= 0) { continue; }
        const Size index = peak - experimental.begin();
        // A decharged graph node has one assigned charge. Unknown peaks alone
        // receive the same bounded hypotheses used during route generation.
        if (! prepared.backgrounds[index].contains(group.charge)) { continue; }
        const double weighting_scale = options.mass_accuracy_tolerance > 0 ? options.mass_accuracy_tolerance : options.tolerance;
        const double weighting_delta = weighting_scale * (options.tolerance_ppm ? ion.mz * 1e-6 : 1.0);
        const double mass_error = peak->getMZ() - ion.mz;
        const double error = weighting_delta > 0 ? mass_error / weighting_delta : (mass_error == 0 ? 0 : std::numeric_limits<double>::infinity());
        by_series[group.series].push_back(
          {Size(peak - experimental.begin()), ion.position, group.charge, peak->getIntensity() * std::exp(-2 * error * error),
           (peak->getMZ() - carrier) * group.charge,
           (group.series == 'A' ? "a" : std::string(1, group.series)) + std::to_string(ion.position) + (group.series == 'A' ? "-B" : "")});
      }
    }
  }
  for (auto& [series, matches] : by_series)
  {
    std::sort(matches.begin(), matches.end(),
              [](const auto& a, const auto& b) { return std::tie(a.position, a.peak, a.charge) < std::tie(b.position, b.peak, b.charge); });
    // Peptide annotation and RNA routes use the same ordinal-series statistic.
    // Ignoring anchoring and peak reuse gives an upper bound: this can only
    // skip series that cannot meet the requested number of route nodes.
    if (matches.size() < options.min_nodes
        || SpectrumAnnotator::computeLongestIonSeries(matches.begin(), matches.end(), options.max_step - 1, [](const auto& match) {
             return match.position;
           }) < options.min_nodes)
    {
      continue;
    }
    std::vector<Route> paths(matches.size());
    Route selected;
    for (Size i = 0; i < matches.size(); ++i)
    {
      const auto& current = matches[i];
      if (current.position <= 2) { paths[i] = {{current}, current.position, current.weighted}; }
      for (Size j = 0; j < i; ++j)
      {
        if (paths[j].matches.empty() || matches[j].position >= current.position || current.position - matches[j].position > options.max_step
            || current.neutral <= matches[j].neutral
            || std::any_of(paths[j].matches.begin(), paths[j].matches.end(), [&](const auto& match) { return match.peak == current.peak; }))
        {
          continue;
        }
        if (paths[j].matches.size() + 1 > paths[i].matches.size()
            || (paths[j].matches.size() + 1 == paths[i].matches.size()
                && paths[j].weighted_intensity + current.weighted > paths[i].weighted_intensity))
        {
          paths[i] = paths[j];
          paths[i].matches.push_back(current);
          paths[i].length = current.position;
          paths[i].weighted_intensity += current.weighted;
        }
      }
      if (paths[i].matches.size() >= options.min_nodes)
      {
        if (paths[i].matches.size() > selected.matches.size()
            || (paths[i].matches.size() == selected.matches.size() && paths[i].weighted_intensity > selected.weighted_intensity))
        {
          selected = paths[i];
        }
      }
    }
    if (! selected.matches.empty())
    {
      selected.series = series;
      selected.scaffold = scaffold;
      result.routes.push_back(std::move(selected));
    }
  }
  return result;
}

RNARouteScorer::Result RNARouteScorer::score(const std::vector<const PreparedRoutes*>& evidence,
                                             const MSSpectrum& experimental,
                                             const PreparedSpectrum& prepared,
                                             bool negative_mode)
{
  Result best;
  if (! experimental.isSorted() || prepared.backgrounds.size() != experimental.size()) { return best; }
  std::map<char, const Route*> selected;
  for (const auto* source : evidence)
  {
    if (! source) { continue; }
    for (const auto& route : source->routes)
    {
      const auto found = selected.find(route.series);
      if (found == selected.end() || route.matches.size() > found->second->matches.size()
          || (route.matches.size() == found->second->matches.size() && route.weighted_intensity > found->second->weighted_intensity))
      {
        selected.insert_or_assign(route.series, &route);
      }
    }
  }
  const Route empty;
  std::vector<const Route*> prefixes {&empty}, suffixes {&empty};
  for (const auto& [series, route] : selected)
  {
    (series >= 'w' ? suffixes : prefixes).push_back(route);
  }
  // Local background uses the same neutral-mass charge hypotheses as the route.
  // Union by original observed peak prevents charge/isobar inflation of intensity.
  for (const auto& prefix : prefixes)
  {
    for (const auto& suffix : suffixes)
    {
      std::map<Size, MatchedNode> unique;
      for (const auto* route : {prefix, suffix})
      {
        for (const auto& match : route->matches)
        {
          const auto found = unique.find(match.peak);
          if (found == unique.end() || match.weighted > found->second.weighted) { unique.insert_or_assign(match.peak, match); }
        }
      }
      if (unique.empty()) { continue; }
      Result result;
      result.nodes = unique.size();
      result.prefix_length = prefix->length;
      result.prefix_scaffold = prefix->scaffold;
      result.suffix_length = suffix->length;
      result.suffix_scaffold = suffix->scaffold;
      std::set<Size> background;
      for (const auto& [index, match] : unique)
      {
        result.matched_intensity += match.weighted;
        const auto& local = prepared.backgrounds[index].at(match.charge);
        background.insert(local.begin(), local.end());
        background.insert(index); // Include evidence even when outside the local top-N.
        PeptideHit::PeakAnnotation annotation;
        annotation.mz = experimental[index].getMZ();
        annotation.intensity = experimental[index].getIntensity();
        annotation.charge = negative_mode ? -Int(match.charge) : Int(match.charge);
        annotation.annotation = match.annotation;
        result.annotations.push_back(std::move(annotation));
      }
      for (Size index : background)
      {
        result.local_intensity += experimental[index].getIntensity();
      }
      result.score = result.nodes * (1 + result.matched_intensity / result.local_intensity);
      if (result.score > best.score) { best = std::move(result); }
    }
  }
  // Equivalent phosphate allocations can produce the same fragment masses via
  // different arithmetic paths. Keep such hypotheses tied instead of selecting
  // a terminal isomer on sub-precision floating-point differences.
  best.score = std::round(best.score * 1e8) / 1e8;
  return best;
}
} // namespace OpenMS

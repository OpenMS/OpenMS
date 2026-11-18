// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#include "RNARouteCandidateGenerator.h"

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string_view>
#include <tuple>
#include <unordered_set>

namespace OpenMS
{
namespace
{
  Size baseIndex(const Ribonucleotide* residue)
  {
    const char base = residue->isModified() ? residue->getOrigin() : residue->getCode()[0];
    return std::string_view("ACGU").find(base);
  }
} // namespace

RNARouteCandidateGenerator::RNARouteCandidateGenerator(const std::vector<NASequence>& database,
                                                       const std::set<const Ribonucleotide*>& modifications,
                                                       const Options& options,
                                                       const std::vector<Byte>& terminal_states):
    database_(database),
    modifications_(modifications),
    terminal_states_(terminal_states),
    options_(options)
{
  if (! std::isfinite(options.tolerance) || options.tolerance <= 0 || ! std::isfinite(options.precursor_tolerance) || options.precursor_tolerance <= 0
      || options.max_length == 0 || options.max_length > 32 || options.max_states == 0 || options.max_completion_nodes == 0
      || options.max_candidates == 0 || options.max_supports == 0 || options.min_route_nodes == 0)
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid RNA route candidate options");
  }
  if (! terminal_states.empty() && terminal_states.size() != database.size())
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Route terminal-state masks must match the database size");
  }
  const double backbone = EmpiricalFormula("H-1PO2").getMonoWeight();
  // The mass graph uses constant residue increments. It omits linkage
  // modifications, whose mass depends on the neighbor: reject them here
  // instead of failing on the first spectrum.
  const auto rejectLinkage = [](const Ribonucleotide& residue) {
    if (residue.getCode().empty() || residue.getCode().back() == '*')
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "graph_routes does not support linkage modifications (codes ending in '*'): '" + residue.getCode() + "'");
    }
  };
  for (const auto* residue : modifications)
  {
    if (! residue || residue->getTermSpecificity() != Ribonucleotide::ANYWHERE)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "graph_routes requires internal variable modifications");
    }
    rejectLinkage(*residue);
  }
  for (const auto& sequence : database)
  {
    for (const auto& residue : sequence)
    {
      rejectLinkage(residue);
    }
    if ((sequence.hasFivePrimeMod() && sequence.getFivePrimeMod()->getCode() != "5'-p")
        || (sequence.hasThreePrimeMod() && sequence.getThreePrimeMod()->getCode() != "3'-p" && sequence.getThreePrimeMod()->getCode() != "3'-c"))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "graph_routes supports OH, phosphate and 3' cyclic phosphate termini");
    }
    double minimum = -backbone;
    for (const auto& residue : sequence)
    {
      double mass = residue.getFormula().getMonoWeight();
      if (! residue.isModified())
      {
        for (const auto* alternative : modifications)
        {
          if (alternative->getTermSpecificity() == Ribonucleotide::ANYWHERE && alternative->getOrigin() == residue.getCode()[0])
          {
            mass = std::min(mass, alternative->getFormula().getMonoWeight());
          }
        }
      }
      minimum += mass + backbone;
    }
    minimum_masses_.push_back(minimum);
  }
  for (Size side = 0; side < 2; ++side)
  {
    auto& trie = tries_[side];
    trie.emplace_back();
    for (Size id = 0; id < database.size(); ++id)
    {
      Size node = 0;
      for (Size pos = 0; pos < database[id].size(); ++pos)
      {
        const auto* residue = database[id][side ? database[id].size() - pos - 1 : pos];
        const Size base = baseIndex(residue);
        if (base >= 4)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Route database requires A/C/G/U residue origins");
        }
        if (trie[node].next[base] < 0)
        {
          const Size next = trie.size();
          trie.emplace_back();
          trie[node].next[base] = next;
        }
        node = trie[node].next[base];
        trie[node].descendants.push_back(id);
      }
    }
  }
}

RNARouteCandidateGenerator::Result RNARouteCandidateGenerator::generate(const MSSpectrum& spectrum,
                                                                        const std::vector<double>& precursor_masses,
                                                                        const RNAModificationDiagnostics::Evidence& evidence,
                                                                        bool negative_mode,
                                                                        Size max_charge,
                                                                        const std::atomic<bool>* cancelled) const
{
  Result result;
  const auto fail = [&](const std::string& reason) {
    result.complete = false;
    result.reason = reason;
    result.candidates.clear();
  };
  const auto cancelledNow = [&]() {
    if (! cancelled || ! cancelled->load(std::memory_order_relaxed)) { return false; }
    fail("cancelled_after_incomplete_search");
    return true;
  };
  if (cancelledNow()) { return result; }
  if (precursor_masses.empty() || spectrum.empty()) { return result; }
  std::vector<std::pair<double, double>> mass_windows;
  const double fraction = options_.precursor_tolerance * 1e-6;
  for (double precursor : precursor_masses)
  {
    if (! std::isfinite(precursor))
    {
      fail("invalid_precursor_mass");
      return result;
    }
    double low = precursor - options_.precursor_tolerance, high = precursor + options_.precursor_tolerance;
    if (options_.precursor_tolerance_ppm)
    {
      // The exact test is |P-M| <= tolerance*M, not tolerance*P.
      // Invert it before querying theoretical candidate-mass intervals.
      low = fraction < 1 ? precursor / (1 + fraction) : -std::numeric_limits<double>::infinity();
      high = fraction < 1 ? precursor / (1 - fraction) : std::numeric_limits<double>::infinity();
    }
    mass_windows.emplace_back(low - 1e-8, high + 1e-8);
  }
  std::sort(mass_windows.begin(), mass_windows.end());
  std::vector<std::pair<double, double>> merged_windows;
  for (const auto& window : mass_windows)
  {
    if (! merged_windows.empty() && window.first <= merged_windows.back().second)
    {
      merged_windows.back().second = std::max(merged_windows.back().second, window.second);
    }
    else
    {
      merged_windows.push_back(window);
    }
  }
  const auto massCompatible = [&](double minimum, double maximum) {
    const auto found = std::lower_bound(merged_windows.begin(), merged_windows.end(), minimum - 1e-8,
                                        [](const auto& window, double value) { return window.second < value; });
    return found != merged_windows.end() && found->first <= maximum + 1e-8;
  };
  if (max_charge == 0)
  {
    fail("invalid_fragment_charge");
    return result;
  }
  auto* db = RibonucleotideDB::getInstance();
  std::set<const Ribonucleotide*> alphabet;
  std::set<const Ribonucleotide*> permitted_variable;
  for (const auto& code : {"A", "C", "G", "U"})
  {
    alphabet.insert(db->getRibonucleotide(code));
  }
  for (const auto* residue : modifications_)
  {
    if (residue->getTermSpecificity() != Ribonucleotide::ANYWHERE)
    {
      fail("terminal_variable_modification");
      return result;
    }
    if (! evidence.tested.contains(residue->getCode()) || evidence.allowed.contains(residue->getCode()))
    {
      alphabet.insert(residue);
      permitted_variable.insert(residue);
    }
  }
  // Fixed/premodified residues have explicit identities and remain available.
  for (const auto& sequence : database_)
  {
    for (const auto& residue : sequence)
    {
      if (residue.isModified()) { alphabet.insert(&residue); }
    }
  }
  const auto original = RNASpectrumGraph::selectPeakIndices(spectrum, options_.max_peaks);
  const double carrier = negative_mode ? -Constants::PROTON_MASS_U : Constants::PROTON_MASS_U;
  const double upper_mass = *std::max_element(precursor_masses.begin(), precursor_masses.end()) + 100;
  MSSpectrum neutral;
  MSSpectrum::IntegerDataArray peak_ids;
  peak_ids.setName("NASE:original_peak");
  MSSpectrum::IntegerDataArray peak_charges;
  peak_charges.setName("NASE:neutral_peak_charge");
  const auto* legacy_charges = RNASpectrumGraph::legacyCharges(spectrum);
  // Unassigned legacy peaks retain bounded charge hypotheses for graph use only.
  // The scoring spectrum is never decharged or mutated by this operation.
  for (Size index : original)
  {
    if (! std::isfinite(spectrum[index].getMZ()) || ! std::isfinite(spectrum[index].getIntensity()) || spectrum[index].getIntensity() <= 0)
    {
      continue;
    }
    std::set<Size> hypotheses;
    const Size assigned = legacy_charges ? Size(std::abs(Int64((*legacy_charges)[index]))) : 0;
    for (Size charge = 1; ! assigned && charge <= max_charge; ++charge)
    {
      hypotheses.insert(charge);
    }
    // Reuse assigned charges exclusively, including those above the bound for
    // unknown peaks. Reinterpreting assigned peaks creates spurious ladders.
    if (assigned) { hypotheses.insert(assigned); }
    for (Size charge : hypotheses)
    {
      const double mass = (spectrum[index].getMZ() - carrier) * charge;
      if (mass > 0 && mass < upper_mass)
      {
        neutral.push_back(Peak1D(mass + Constants::PROTON_MASS_U, spectrum[index].getIntensity()));
        peak_ids.push_back(index);
        peak_charges.push_back(charge);
      }
    }
  }
  neutral.getIntegerDataArrays().push_back(std::move(peak_ids));
  neutral.getIntegerDataArrays().push_back(std::move(peak_charges));
  neutral.sortByPosition();
  std::vector<double> neutral_tolerances;
  double seed_tolerance = options_.tolerance;
  if (! options_.tolerance_ppm)
  {
    // The tolerance is specified in observed m/z. Decharging scales each
    // peak's absolute uncertainty by its own charge, including graph edges.
    for (Int charge : neutral.getIntegerDataArrays()[1])
    {
      neutral_tolerances.push_back(options_.tolerance * charge);
      seed_tolerance = std::max(seed_tolerance, neutral_tolerances.back());
    }
  }
  RNASpectrumGraph::Options graph_options;
  graph_options.tolerance = options_.tolerance;
  graph_options.tolerance_ppm = options_.tolerance_ppm;
  graph_options.min_tag_length = 1;
  graph_options.max_tag_length = options_.max_length;
  graph_options.max_peaks = std::max<Size>(2, neutral.size());
  graph_options.max_charge = 1;
  graph_options.max_edges = options_.max_states;
  const auto graph = RNASpectrumGraph(alphabet, graph_options).build(neutral, 1, neutral_tolerances);
  result.graph_nodes = graph.peaks;
  result.graph_edges = graph.edges;
  if (! graph.fallback_reason.empty())
  {
    if (graph.fallback_reason != "no_tags") { fail("graph:" + graph.fallback_reason); }
    return result;
  }
  std::vector<const Ribonucleotide*> residues(graph.alphabet->code_ids.size());
  for (const auto* residue : alphabet)
  {
    residues[graph.alphabet->code_ids.at(residue->getCode())] = residue;
  }
  struct Route
  {
    Size side, terminal, trie_node;
    std::vector<UInt> codes;
    double delta_mass = 0;
  };
  std::vector<Route> routes;
  std::set<std::tuple<Size, Size, std::vector<UInt>>> route_keys;
  std::set<std::tuple<Size, Size, Size, Size, std::vector<UInt>, std::vector<Size>, Size>> states;
  std::vector<Size> last_peak_node(spectrum.size(), 0);
  for (Size node = 0; node < graph.peak_indices.size(); ++node)
  {
    last_peak_node[neutral.getIntegerDataArrays()[0][graph.peak_indices[node]]] = node;
  }
  const double backbone = EmpiricalFormula("H-1PO2").getMonoWeight();
  const double phosphate = EmpiricalFormula("HPO3").getMonoWeight();
  const double water = EmpiricalFormula("H2O").getMonoWeight();
  std::vector<double> residue_deltas(residues.size());
  for (Size i = 0; i < residues.size(); ++i)
  {
    const Size base = baseIndex(residues[i]);
    if (base < 4)
    {
      residue_deltas[i] = residues[i]->getFormula().getMonoWeight()
                          - db->getRibonucleotide(std::string(1, std::string_view("ACGU")[base]))->getFormula().getMonoWeight();
    }
  }
  const std::array<double, 3> terminal_mass {0, phosphate, phosphate - water};
  const NASequence reference = NASequence::fromString("A");
  const double reference_increment = db->getRibonucleotide("A")->getFormula().getMonoWeight() + backbone;
  const auto advance = [&](Size node, const RNASpectrumGraph::Transition& step, Size side) {
    for (Size i = 0; i < step.length; ++i)
    {
      const Size base = baseIndex(residues[step.residues[i]]);
      if (base >= 4 || tries_[side][node].next[base] < 0) { return Size(-1); }
      node = tries_[side][node].next[base];
    }
    return node;
  };
  const auto within = [&](Size peak, double expected) {
    const double tolerance = options_.tolerance_ppm ? options_.tolerance * expected * 1e-6 : neutral_tolerances[peak];
    return std::abs(neutral[peak].getMZ() - expected) <= tolerance;
  };
  // Equivalent mass offsets share state, but terminal assignments stay distinct.
  for (Size side = 0; side < 2 && result.complete; ++side)
  {
    for (Size terminal = 0; terminal < (side ? 3 : 2) && result.complete; ++terminal)
    {
      for (auto type : {NASequence::AIon, NASequence::BIon, NASequence::CIon, NASequence::DIon})
      {
        const double offset = reference.getMonoWeight(type, 0) - reference_increment + terminal_mass[terminal] + Constants::PROTON_MASS_U;
        // Prefix a/b/c/d and suffix z/y/x/w have the same OH reference offsets.
        std::function<void(Size, Size, std::vector<UInt>&, double, Size, std::vector<Size>)> extend;
        extend = [&](Size node, Size trie_node, std::vector<UInt>& codes, double increment, Size observed_nodes, std::vector<Size> used_peaks) {
          if (! result.complete || cancelledNow()) { return; }
          const Size original_peak = neutral.getIntegerDataArrays()[0][graph.peak_indices[node]];
          const auto position = std::lower_bound(used_peaks.begin(), used_peaks.end(), original_peak);
          if (position != used_peaks.end() && *position == original_peak) { return; }
          used_peaks.insert(position, original_peak);
          if (options_.compress_route_states)
          {
            // Once all charge hypotheses for a used peak lie behind this DAG
            // node, its identity cannot constrain a future extension. Forget
            // only those expired identities, preserving exact future choices.
            std::erase_if(used_peaks, [&](Size peak) { return last_peak_node[peak] <= node; });
          }
          // Future compatibility depends on which original peaks are already
          // used. Distinct charge hypotheses of one peak cannot supply two rungs.
          const auto key = std::tuple(side, terminal, Size(type), node, codes, used_peaks, std::min(observed_nodes, options_.min_route_nodes));
          if (! states.insert(key).second) { return; }
          if (++result.route_states > options_.max_states)
          {
            fail("route_state_limit");
            return;
          }
          if (observed_nodes >= options_.min_route_nodes && route_keys.emplace(side, terminal, codes).second)
          {
            double delta = 0;
            for (UInt code : codes)
            {
              delta += residue_deltas[code];
            }
            routes.push_back({side, terminal, trie_node, codes, delta});
          }
          for (const auto& edge : graph.layers[0].outgoing[node])
          {
            if (codes.size() + edge.length > options_.max_length) { continue; }
            for (Size step_id = edge.first_transition; step_id < edge.past_transition; ++step_id)
            {
              const auto& step = graph.alphabet->transitions[step_id];
              if (step.length != edge.length) { continue; }
              const Size next_trie = advance(trie_node, step, side);
              if (next_trie == Size(-1)) { continue; }
              const double next_mass = increment + step.mass;
              if (! within(graph.peak_indices[edge.to], next_mass + offset)) { continue; }
              const Size old_size = codes.size();
              codes.insert(codes.end(), step.residues.begin(), step.residues.begin() + step.length);
              extend(edge.to, next_trie, codes, next_mass, observed_nodes + 1, used_peaks);
              codes.resize(old_size);
              if (! result.complete || cancelledNow()) { return; }
            }
          }
        };
        for (const auto& seed : graph.alphabet->transitions)
        {
          const Size trie_node = advance(0, seed, side);
          if (trie_node == Size(-1)) { continue; }
          const double expected = seed.mass + offset;
          const double delta = options_.tolerance_ppm ? options_.tolerance * expected * 1e-6 : seed_tolerance;
          auto first = std::lower_bound(graph.peak_indices.begin(), graph.peak_indices.end(), expected - delta,
                                        [&](Size index, double mz) { return neutral[index].getMZ() < mz; });
          for (auto it = first; it != graph.peak_indices.end() && neutral[*it].getMZ() <= expected + delta; ++it)
          {
            if (! within(*it, expected)) { continue; }
            std::vector<UInt> codes(seed.residues.begin(), seed.residues.begin() + seed.length);
            extend(it - graph.peak_indices.begin(), trie_node, codes, seed.mass, 1, {});
          }
          if (! result.complete) { break; }
        }
      }
    }
  }
  result.routes = routes.size();
  if (! result.complete || routes.empty()) { return result; }
  std::vector<std::array<std::vector<Size>, 2>> by_database(database_.size());
  for (Size i = 0; i < routes.size(); ++i)
  {
    const auto& route = routes[i];
    for (Size id : tries_[route.side][route.trie_node].descendants)
    {
      if (massCompatible(minimum_masses_[id], std::numeric_limits<double>::infinity())) { by_database[id][route.side].push_back(i); }
    }
  }
  const std::array<std::pair<Size, Size>, 6> terminal_states {{{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}}};
  std::map<std::tuple<Size, std::string, Size>, Size> emitted;
  std::map<std::string, std::shared_ptr<const Scaffold>> scaffold_cache;
  Size support_states = 0;
  constexpr double mass_quantum = 1e-5;
  constexpr Size max_mass_bins = 4096, max_cached_bins = 1000000;
  std::map<std::vector<Int64>, Size> mass_pools;
  std::map<std::vector<Size>, std::vector<Int64>> attainable_cache {{{}, {0}}};
  Size cached_bins = 1;
  for (Size id = 0; id < database_.size() && result.complete; ++id)
  {
    if (cancelledNow()) { return result; }
    auto& supported = by_database[id];
    if (supported[0].empty() && supported[1].empty()) { continue; }
    if (options_.compress_equivalent_routes)
    {
      // A longer route carrying the same placed modifications only adds
      // canonical constraints. With either labile or generic completion, its shorter
      // witnessed subroute already permits every completion. Keep the least
      // restrictive witness for each modification mask and terminal state.
      for (Size side = 0; side < 2; ++side)
      {
        std::map<std::pair<Size, std::vector<std::pair<Size, UInt>>>, Size> masks;
        for (Size route_id : supported[side])
        {
          const auto& route = routes[route_id];
          std::vector<std::pair<Size, UInt>> modifications;
          for (Size pos = 0; pos < route.codes.size(); ++pos)
          {
            if (residues[route.codes[pos]]->isModified()) { modifications.emplace_back(pos, route.codes[pos]); }
          }
          const auto key = std::pair(route.terminal, modifications);
          const auto found = masks.find(key);
          if (found == masks.end() || route.codes.size() < routes[found->second].codes.size()) { masks.insert_or_assign(key, route_id); }
        }
        supported[side].clear();
        for (const auto& [key, route_id] : masks)
        {
          supported[side].push_back(route_id);
        }
      }
    }
    supported[0].push_back(Size(-1));
    supported[1].push_back(Size(-1));
    const Size length = database_[id].size();
    const bool indexed = options_.index_precursor_mass && length <= 32
                         && std::none_of(database_[id].begin(), database_[id].end(), [](const auto& residue) { return residue.isModified(); });
    NASequence bare = database_[id];
    bare.setFivePrimeMod(nullptr);
    bare.setThreePrimeMod(nullptr);
    const double bare_mass = bare.getMonoWeight();
    double completion_delta_min = 0, completion_delta_max = 0;
    std::vector<double> completion_minima, completion_maxima;
    for (const auto& original_residue : database_[id])
    {
      double smallest = 0, largest = 0;
      for (const auto* residue : permitted_variable)
      {
        if (residue->getOrigin() == original_residue.getCode()[0]
            && (options_.complete_unanchored_modifications || evidence.labile.contains(residue->getCode())))
        {
          const double delta = residue->getFormula().getMonoWeight() - original_residue.getFormula().getMonoWeight();
          smallest = std::min(smallest, delta);
          largest = std::max(largest, delta);
        }
      }
      // Include even occupied positions: this broadens the bound safely.
      completion_minima.push_back(smallest);
      completion_maxima.push_back(largest);
    }
    // Only unanchored completions consume this budget. Intact route deltas
    // remain exact, and occupied-site deltas merely widen these intervals.
    std::sort(completion_minima.begin(), completion_minima.end());
    std::sort(completion_maxima.begin(), completion_maxima.end(), std::greater<double>());
    const Size completion_budget = std::min(length, options_.max_unanchored_modifications);
    for (Size pos = 0; pos < completion_budget; ++pos)
    {
      completion_delta_min += completion_minima[pos];
      completion_delta_max += completion_maxima[pos];
    }
    const auto routeLength = [&](Size route) { return route == Size(-1) ? 0 : routes[route].codes.size(); };
    const auto routeDelta = [&](Size route) { return route == Size(-1) ? 0.0 : routes[route].delta_mass; };
    std::vector<Size> sorted_suffixes = supported[1];
    std::sort(sorted_suffixes.begin(), sorted_suffixes.end(),
              [&](Size a, Size b) { return std::pair(routeDelta(a), a) < std::pair(routeDelta(b), b); });
    std::vector<std::vector<Size>> suffixes_by_length(length + 1);
    for (Size route : supported[1])
    {
      if (routeLength(route) <= length) { suffixes_by_length[routeLength(route)].push_back(route); }
    }
    std::map<std::tuple<Size, Size, std::vector<const Ribonucleotide*>>, std::set<std::vector<Size>>> patterns;
    for (Size prefix_id : supported[0])
    {
      std::vector<std::pair<Size, Size>> pairs;
      for (Size state = 0; state < terminal_states.size(); ++state)
      {
        const auto [five, three] = terminal_states[state];
        if (! terminal_states_.empty() && ! ((terminal_states_[id] >> state) & 1)) { continue; }
        if (prefix_id != Size(-1) && routes[prefix_id].terminal != five) { continue; }
        if (! indexed)
        {
          for (Size suffix_id : supported[1])
          {
            pairs.emplace_back(suffix_id, state);
          }
          continue;
        }
        const Size prefix_length = routeLength(prefix_id);
        if (prefix_length > length) { continue; }
        const double origin = bare_mass + terminal_mass[five] + terminal_mass[three] + routeDelta(prefix_id);
        std::vector<std::pair<double, double>> intervals;
        for (const auto& window : merged_windows)
        {
          intervals.emplace_back(window.first - origin - completion_delta_max - 1e-8, window.second - origin - completion_delta_min + 1e-8);
        }
        // Subtracting the same labile bounds preserves order, but can merge windows.
        std::vector<std::pair<double, double>> merged;
        for (const auto& interval : intervals)
        {
          if (! merged.empty() && interval.first <= merged.back().second) { merged.back().second = std::max(merged.back().second, interval.second); }
          else
          {
            merged.push_back(interval);
          }
        }
        for (const auto& interval : merged)
        {
          auto first = std::lower_bound(sorted_suffixes.begin(), sorted_suffixes.end(), interval.first,
                                        [&](Size route, double value) { return routeDelta(route) < value; });
          for (auto it = first; it != sorted_suffixes.end() && routeDelta(*it) <= interval.second; ++it)
          {
            if (prefix_length + routeLength(*it) <= length) { pairs.emplace_back(*it, state); }
          }
        }
        // Overlapping routes share placements, so their deltas are not additive.
        // Retain the exhaustive compatibility check for these pairs.
        // Full-length overlaps can add distinct labile route witnesses even
        // when a one-sided route already permits the same completed sequence.
        for (Size suffix_length = length - prefix_length + 1; suffix_length <= length; ++suffix_length)
        {
          for (Size suffix_id : suffixes_by_length[suffix_length])
          {
            pairs.emplace_back(suffix_id, state);
          }
        }
      }
      for (const auto& [suffix_id, state] : pairs)
      {
        if (prefix_id == Size(-1) && suffix_id == Size(-1)) { continue; }
        const auto [five, three] = terminal_states[state];
        // A mask already accounts for the stored termini and where the enzyme replaces them.
        if (terminal_states_.empty()
            && ((database_[id].hasFivePrimeMod() && five != 1)
                || (database_[id].hasThreePrimeMod() && three != (database_[id].getThreePrimeMod()->getCode() == "3'-p" ? 1 : 2))))
        {
          continue;
        }
        if ((prefix_id != Size(-1) && routes[prefix_id].terminal != five) || (suffix_id != Size(-1) && routes[suffix_id].terminal != three))
        {
          continue;
        }
        std::vector<const Ribonucleotide*> fixed(database_[id].size(), nullptr);
        bool compatible = true;
        for (Size side = 0; side < 2 && compatible; ++side)
        {
          const Size route_id = side ? suffix_id : prefix_id;
          if (route_id == Size(-1)) { continue; }
          const auto& route = routes[route_id];
          if (route.codes.size() > fixed.size() || (! options_.allow_full_length_routes && route.codes.size() == fixed.size()))
          {
            compatible = false;
            break;
          }
          for (Size pos = 0; pos < route.codes.size(); ++pos)
          {
            const Size location = side ? fixed.size() - pos - 1 : pos;
            const auto* residue = residues[route.codes[pos]];
            if (database_[id][location]->isModified() && database_[id][location] != residue)
            {
              compatible = false;
              break;
            }
            // Premodified residues belong to their original database sites.
            // Moving one to another site requires an explicitly permitted variable mod.
            if (residue->isModified() && database_[id][location] != residue && ! permitted_variable.contains(residue))
            {
              compatible = false;
              break;
            }
            if (fixed[location] && fixed[location] != residue)
            {
              const auto* modified = residue->isModified() ? residue : fixed[location];
              const auto* canonical = residue->isModified() ? fixed[location] : residue;
              // A canonical fragment route may coexist with an intact labile
              // route at the same site. The full precursor retains the labile
              // residue; each side keeps its own witnessed fragment scaffold.
              if (canonical->isModified() || ! evidence.labile.contains(modified->getCode()) || baseIndex(canonical) != baseIndex(modified))
              {
                compatible = false;
                break;
              }
              fixed[location] = modified;
            }
            else
            {
              fixed[location] = residue;
            }
          }
        }
        if (! compatible) { continue; }
        if (! options_.complete_unanchored_modifications && options_.compress_equivalent_routes)
        {
          // An uncovered canonical site has exactly the same canonical plus
          // labile choices as an explicitly canonical route site.
          for (auto& residue : fixed)
          {
            if (residue && ! residue->isModified()) { residue = nullptr; }
          }
        }
        auto& masks = patterns[std::tuple(five, three, fixed)];
        for (Size side = 0; side < 2; ++side)
        {
          const Size route_id = side ? suffix_id : prefix_id;
          if (route_id == Size(-1)) { continue; }
          std::vector<Size> intact;
          for (Size pos = 0; pos < routes[route_id].codes.size(); ++pos)
          {
            if (evidence.labile.contains(residues[routes[route_id].codes[pos]]->getCode())) { intact.push_back(side ? fixed.size() - pos - 1 : pos); }
          }
          std::sort(intact.begin(), intact.end());
          if (masks.insert(std::move(intact)).second && ++support_states > options_.max_supports)
          {
            fail("support_limit");
            return result;
          }
        }
      }
    }
    // Complete each unique placement pattern once after collecting all of its
    // side-specific labile witnesses. A duplicate pattern can add evidence.
    for (const auto& pattern_entry : patterns)
    {
      if (cancelledNow()) { return result; }
      const auto& support_masks = pattern_entry.second;
      const auto five = std::get<0>(pattern_entry.first);
      const auto three = std::get<1>(pattern_entry.first);
      const auto& fixed = std::get<2>(pattern_entry.first);
      std::vector<std::vector<const Ribonucleotide*>> choices(fixed.size());
      for (Size pos = 0; pos < fixed.size(); ++pos)
      {
        const auto* original_residue = database_[id][pos];
        choices[pos].push_back(fixed[pos] ? fixed[pos] : original_residue);
        if (original_residue->isModified()) { continue; }
        for (const auto* residue : permitted_variable)
        {
          if (! residue->isModified() || residue->getOrigin() != original_residue->getCode()[0]) { continue; }
          if ((! fixed[pos] && options_.complete_unanchored_modifications)
              || ((! fixed[pos] || ! fixed[pos]->isModified()) && evidence.labile.contains(residue->getCode())))
          {
            choices[pos].push_back(residue);
          }
        }
        std::sort(choices[pos].begin(), choices[pos].end());
        choices[pos].erase(std::unique(choices[pos].begin(), choices[pos].end()), choices[pos].end());
      }
      std::vector<double> low(choices.size() + 1, 0), high(choices.size() + 1, 0);
      std::vector<std::vector<double>> mono_masses(choices.size());
      for (Size pos = choices.size(); pos-- > 0;)
      {
        double smallest = std::numeric_limits<double>::infinity(), largest = -smallest;
        for (const auto* residue : choices[pos])
        {
          const double mono_mass = residue->getFormula().getMonoWeight();
          mono_masses[pos].push_back(mono_mass);
          const double mass = mono_mass + backbone;
          smallest = std::min(smallest, mass);
          largest = std::max(largest, mass);
        }
        low[pos] = low[pos + 1] + smallest;
        high[pos] = high[pos + 1] + largest;
      }
      const auto unanchoredCost = [&](Size pos, const Ribonucleotide* residue) -> Size {
        return (! fixed[pos] || ! fixed[pos]->isModified()) && residue->isModified() && ! database_[id][pos]->isModified();
      };
      std::vector<std::vector<double>> budget_low, budget_high;
      if (options_.index_precursor_mass && options_.max_unanchored_modifications < choices.size())
      {
        const Size budget = options_.max_unanchored_modifications;
        budget_low.assign(choices.size() + 1, std::vector<double>(budget + 1, 0));
        budget_high = budget_low;
        for (Size pos = choices.size(); pos-- > 0;)
        {
          for (Size remaining = 0; remaining <= budget; ++remaining)
          {
            double smallest = std::numeric_limits<double>::infinity(), largest = -smallest;
            for (Size choice = 0; choice < choices[pos].size(); ++choice)
            {
              const Size cost = unanchoredCost(pos, choices[pos][choice]);
              if (cost > remaining) { continue; }
              const double mass = mono_masses[pos][choice] + backbone;
              smallest = std::min(smallest, mass + budget_low[pos + 1][remaining - cost]);
              largest = std::max(largest, mass + budget_high[pos + 1][remaining - cost]);
            }
            budget_low[pos][remaining] = smallest;
            budget_high[pos][remaining] = largest;
          }
        }
      }
      const auto lowerSuffix
        = [&](Size pos, Size used) { return budget_low.empty() ? low[pos] : budget_low[pos][options_.max_unanchored_modifications - used]; };
      const auto upperSuffix
        = [&](Size pos, Size used) { return budget_high.empty() ? high[pos] : budget_high[pos][options_.max_unanchored_modifications - used]; };
      std::vector<const Ribonucleotide*> placed(choices.size());
      const double terminal_offset = -backbone + terminal_mass[five] + terminal_mass[three];
      // Share exact precursor acceptance and scaffold bookkeeping between
      // indexed placement lookup and the exhaustive recursive reference.
      const auto emitCandidate = [&]() {
        NASequence candidate(placed, five ? db->getRibonucleotide("5'-p") : nullptr,
                             three == 1 ? db->getRibonucleotide("3'-p") : (three == 2 ? db->getRibonucleotide("3'-c") : nullptr));
        const double exact = candidate.getMonoWeight();
        const double delta = options_.precursor_tolerance * (options_.precursor_tolerance_ppm ? exact * 1e-6 : 1.0);
        if (std::none_of(precursor_masses.begin(), precursor_masses.end(), [&](double mass) { return std::abs(mass - exact) <= delta; })) { return; }
        const auto label = candidate.toString();
        std::vector<std::shared_ptr<const Scaffold>> scaffold_variants;
        const bool has_labile
          = std::any_of(candidate.cbegin(), candidate.cend(), [&](const auto& residue) { return evidence.labile.contains(residue.getCode()); });
        if (has_labile)
        {
          for (const auto& intact : support_masks)
          {
            auto scaffold_sequence = candidate;
            for (Size location = 0; location < candidate.size(); ++location)
            {
              if (! database_[id][location]->isModified() && evidence.labile.contains(candidate[location]->getCode())
                  && ! std::binary_search(intact.begin(), intact.end(), location))
              {
                scaffold_sequence.set(location, database_[id][location]);
              }
            }
            const auto scaffold_label = scaffold_sequence.toString();
            if (scaffold_label == label) { continue; }
            auto found = scaffold_cache.find(scaffold_label);
            if (found == scaffold_cache.end())
            {
              if (scaffold_cache.size() >= options_.max_supports)
              {
                fail("support_limit");
                return;
              }
              found
                = scaffold_cache.emplace(scaffold_label, std::make_shared<const Scaffold>(Scaffold {std::move(scaffold_sequence), scaffold_label}))
                    .first;
            }
            if (std::find(scaffold_variants.begin(), scaffold_variants.end(), found->second) == scaffold_variants.end())
            {
              scaffold_variants.push_back(found->second);
            }
          }
        }
        for (Size precursor = 0; precursor < precursor_masses.size(); ++precursor)
        {
          if (std::abs(precursor_masses[precursor] - exact) > delta) { continue; }
          const auto key = std::tuple(id, label, precursor);
          auto found = emitted.find(key);
          if (found == emitted.end())
          {
            if (result.candidates.size() >= options_.max_candidates)
            {
              fail("candidate_limit");
              return;
            }
            found = emitted.emplace(key, result.candidates.size()).first;
            result.candidates.push_back({id, precursor, candidate, {}});
          }
          auto& supports = result.candidates[found->second].support_scaffolds;
          for (const auto& scaffold : scaffold_variants)
          {
            if (std::find(supports.begin(), supports.end(), scaffold) != supports.end()) { continue; }
            if (result.support_variants >= options_.max_supports)
            {
              fail("support_limit");
              return;
            }
            supports.push_back(scaffold);
            ++result.support_variants;
          }
        }
      };
      // Each site normally has one zero-cost canonical/fixed choice. With at
      // most two unsupported placements, query exact one-site mass shifts and
      // their complementary pairs instead of traversing every intervening site.
      // The recursive path remains the oracle and handles larger budgets.
      bool direct_lookup = options_.index_precursor_mass && options_.max_unanchored_modifications <= 2;
      std::vector<Size> baseline(choices.size());
      if (direct_lookup)
      {
        for (Size pos = 0; pos < choices.size(); ++pos)
        {
          Size zero_cost = 0;
          for (Size choice = 0; choice < choices[pos].size(); ++choice)
          {
            if (unanchoredCost(pos, choices[pos][choice]) == 0)
            {
              baseline[pos] = choice;
              ++zero_cost;
            }
          }
          if (zero_cost != 1)
          {
            direct_lookup = false;
            break;
          }
        }
      }
      if (direct_lookup)
      {
        const auto tick = [&]() {
          if (! result.complete || cancelledNow()) { return false; }
          if (++result.completion_nodes > options_.max_completion_nodes)
          {
            fail("completion_node_limit");
            return false;
          }
          return true;
        };
        if (! tick()) { return result; }
        if (! massCompatible(lowerSuffix(0, 0) + terminal_offset, upperSuffix(0, 0) + terminal_offset)) { continue; }
        double origin = terminal_offset;
        Size fixed_modifications = 0;
        for (Size pos = 0; pos < choices.size(); ++pos)
        {
          if (! tick()) { return result; }
          placed[pos] = choices[pos][baseline[pos]];
          origin += mono_masses[pos][baseline[pos]] + backbone;
          fixed_modifications += placed[pos]->isModified() && ! database_[id][pos]->isModified();
        }
        if (fixed_modifications > options_.max_modifications) { continue; }
        const Size budget = std::min(options_.max_unanchored_modifications, options_.max_modifications - fixed_modifications);
        struct Placement
        {
          Size position;
          const Ribonucleotide* residue;
          double delta;
        };
        std::vector<Placement> alternatives;
        if (budget)
        {
          for (Size pos = 0; pos < choices.size(); ++pos)
          {
            for (Size choice = 0; choice < choices[pos].size(); ++choice)
            {
              if (choice != baseline[pos])
              {
                alternatives.push_back({pos, choices[pos][choice], mono_masses[pos][choice] - mono_masses[pos][baseline[pos]]});
              }
            }
          }
          std::sort(alternatives.begin(), alternatives.end(), [](const auto& a, const auto& b) {
            return std::tuple(a.delta, a.position, a.residue->getCode()) < std::tuple(b.delta, b.position, b.residue->getCode());
          });
        }
        constexpr double margin = 1e-7;
        if (massCompatible(origin - margin, origin + margin))
        {
          emitCandidate();
          if (! result.complete) { return result; }
        }
        if (! alternatives.empty())
        {
          auto window_it = std::lower_bound(merged_windows.begin(), merged_windows.end(), origin + alternatives.front().delta - margin,
                                            [](const auto& window, double mass) { return window.second < mass; });
          for (; window_it != merged_windows.end() && window_it->first <= origin + alternatives.back().delta + margin; ++window_it)
          {
            if (! tick()) { return result; }
            const auto& window = *window_it;
            const double lower = window.first - origin - margin, upper = window.second - origin + margin;
            auto first = std::lower_bound(alternatives.begin(), alternatives.end(), lower,
                                          [](const auto& placement, double delta) { return placement.delta < delta; });
            for (auto it = first; it != alternatives.end() && it->delta <= upper; ++it)
            {
              if (! tick()) { return result; }
              placed[it->position] = it->residue;
              emitCandidate();
              placed[it->position] = choices[it->position][baseline[it->position]];
              if (! result.complete) { return result; }
            }
          }
          if (budget == 2)
          {
            for (Size group = 0; group + 1 < alternatives.size();)
            {
              const double shift = alternatives[group].delta;
              const auto end = std::upper_bound(alternatives.begin() + group, alternatives.end(), shift,
                                                [](double delta, const auto& placement) { return delta < placement.delta; });
              const Size group_end = end - alternatives.begin();
              std::vector<std::pair<Size, Size>> partner_ranges;
              // Equal shifts share a mass query, while every position and
              // isobaric residue identity still gets its own emitted placement.
              auto paired_window
                = std::lower_bound(merged_windows.begin(), merged_windows.end(), origin + shift + alternatives[group + 1].delta - margin,
                                   [](const auto& window, double mass) { return window.second < mass; });
              for (; paired_window != merged_windows.end() && paired_window->first <= origin + shift + alternatives.back().delta + margin;
                   ++paired_window)
              {
                if (! tick()) { return result; }
                const double lower = paired_window->first - origin - shift - margin, upper = paired_window->second - origin - shift + margin;
                const auto found = std::lower_bound(alternatives.begin(), alternatives.end(), lower,
                                                    [](const auto& placement, double delta) { return placement.delta < delta; });
                const auto last
                  = std::upper_bound(found, alternatives.end(), upper, [](double delta, const auto& placement) { return delta < placement.delta; });
                if (found != last) { partner_ranges.emplace_back(found - alternatives.begin(), last - alternatives.begin()); }
              }
              for (Size i = group; i < group_end; ++i)
              {
                const auto& first = alternatives[i];
                placed[first.position] = first.residue;
                for (const auto& [first_partner, last_partner] : partner_ranges)
                {
                  for (Size j = std::max(i + 1, first_partner); j < last_partner; ++j)
                  {
                    if (! tick()) { return result; }
                    const auto& second = alternatives[j];
                    if (first.position == second.position) { continue; }
                    placed[second.position] = second.residue;
                    emitCandidate();
                    placed[second.position] = choices[second.position][baseline[second.position]];
                    if (! result.complete) { return result; }
                  }
                }
                placed[first.position] = choices[first.position][baseline[first.position]];
              }
              group = group_end;
            }
          }
        }
        continue;
      }
      // A min/max envelope includes many unattainable labile mass sums. Index
      // conservative bins of attainable suffix sums without merging placements.
      // Flooring each variable shift contributes at most one quantum of uncertainty;
      // the widened interval is only a pruning bound, never a mass acceptance.
      // Cache by the multiset of choice pools, independent of fixed residues
      // and placement order. Thus repeated labile-completion patterns share
      // the same small mass lookup while every sequence placement is retained.
      std::vector<const std::vector<Int64>*> attainable(choices.size() + 1, nullptr);
      std::vector<double> reference_mass(choices.size() + 1, 0);
      std::vector<Size> variable_sites(choices.size() + 1, 0);
      if (options_.index_precursor_mass && choices.size() <= 32
          && std::count_if(choices.begin(), choices.end(), [](const auto& values) { return values.size() > 1; }) >= 4
          && massCompatible(lowerSuffix(0, 0) + terminal_offset, upperSuffix(0, 0) + terminal_offset))
      {
        attainable.back() = &attainable_cache.begin()->second;
        std::vector<Size> key;
        for (Size pos = choices.size(); pos-- > 0;)
        {
          reference_mass[pos] = reference_mass[pos + 1] + mono_masses[pos][0] + backbone;
          variable_sites[pos] = variable_sites[pos + 1] + (choices[pos].size() > 1);
          if (choices[pos].size() == 1)
          {
            attainable[pos] = attainable[pos + 1];
            continue;
          }
          std::vector<Int64> shifts;
          bool valid = true;
          for (double mono : mono_masses[pos])
          {
            const double scaled = std::floor((mono - mono_masses[pos][0]) / mass_quantum);
            if (! std::isfinite(scaled) || std::abs(scaled) > double(std::numeric_limits<Int64>::max() / 64))
            {
              valid = false;
              break;
            }
            shifts.push_back(Int64(scaled));
          }
          if (! valid) { break; } // Bounds alone remain safe.
          std::sort(shifts.begin(), shifts.end());
          shifts.erase(std::unique(shifts.begin(), shifts.end()), shifts.end());
          const Size pool = mass_pools.try_emplace(shifts, mass_pools.size()).first->second;
          key.insert(std::upper_bound(key.begin(), key.end(), pool), pool);
          const auto cached = attainable_cache.find(key);
          if (cached != attainable_cache.end())
          {
            if (! cached->second.empty()) { attainable[pos] = &cached->second; }
            continue;
          }
          if (! attainable[pos + 1] || cached_bins >= max_cached_bins || attainable_cache.size() >= 50000) { continue; }
          std::vector<Int64> sums;
          for (Int64 shift : shifts)
          {
            for (Int64 suffix : *attainable[pos + 1])
            {
              sums.push_back(shift + suffix);
            }
          }
          std::sort(sums.begin(), sums.end());
          sums.erase(std::unique(sums.begin(), sums.end()), sums.end());
          if (sums.size() > max_mass_bins || cached_bins + sums.size() > max_cached_bins) { sums.clear(); }
          cached_bins += sums.size();
          const auto stored = attainable_cache.emplace(key, std::move(sums)).first;
          if (! stored->second.empty()) { attainable[pos] = &stored->second; }
        }
      }
      std::function<void(Size, Size, Size, double)> complete;
      complete = [&](Size pos, Size modifications, Size unanchored, double mass) {
        if (! result.complete || cancelledNow()) { return; }
        if (++result.completion_nodes > options_.max_completion_nodes)
        {
          fail("completion_node_limit");
          return;
        }
        if (modifications > options_.max_modifications || unanchored > options_.max_unanchored_modifications) { return; }
        const double minimum = mass + lowerSuffix(pos, unanchored) + terminal_offset, maximum = mass + upperSuffix(pos, unanchored) + terminal_offset;
        if (! massCompatible(minimum, maximum)) { return; }
        if (attainable[pos] && (pos == 0 || (pos < choices.size() && choices[pos].size() > 1)))
        {
          const double width = variable_sites[pos] * mass_quantum + 1e-8;
          // The suffix bins are sorted. Query each precursor interval rather
          // than scanning every attainable sum at every completion node.
          // Widen only this pruning bound to retain boundary-rounding cases;
          // exact precursor validation below remains unchanged.
          const double origin = mass + terminal_offset + reference_mass[pos];
          bool reachable = false;
          for (const auto& window : merged_windows)
          {
            const auto found = std::lower_bound(attainable[pos]->begin(), attainable[pos]->end(), window.first - width - 1e-7,
                                                [&](Int64 sum, double lower) { return origin + sum * mass_quantum < lower; });
            if (found != attainable[pos]->end() && origin + *found * mass_quantum <= window.second + 1e-7)
            {
              reachable = true;
              break;
            }
          }
          if (! reachable) { return; }
        }
        if (pos != choices.size())
        {
          for (Size choice = 0; choice < choices[pos].size(); ++choice)
          {
            const auto* residue = choices[pos][choice];
            const Size cost = unanchoredCost(pos, residue);
            if (cost > options_.max_unanchored_modifications - unanchored) { continue; }
            placed[pos] = residue;
            complete(pos + 1, modifications + (residue->isModified() && ! database_[id][pos]->isModified()), unanchored + cost,
                     mass + mono_masses[pos][choice] + backbone);
            if (! result.complete || cancelledNow()) { return; }
          }
          return;
        }
        emitCandidate();
      };
      complete(0, 0, 0, 0);
      if (! result.complete) { return result; }
    }
  }
  for (auto& candidate : result.candidates)
  {
    std::sort(candidate.support_scaffolds.begin(), candidate.support_scaffolds.end(),
              [](const auto& a, const auto& b) { return a->label < b->label; });
  }
  return result;
}
} // namespace OpenMS

// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "ModifiedNASequenceGenerator.h"

#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

using std::vector;
using std::set;
using std::map;

namespace OpenMS
{
  // static
  void ModifiedNASequenceGenerator::applyFixedModifications(
    const set<ConstRibonucleotidePtr>& fixed_mods,
    NASequence& seq)
  {
    // apply modifications at chain ends
    std::for_each(fixed_mods.begin(), fixed_mods.end(), [&seq] (const ConstRibonucleotidePtr& f)
      {
        if (f->getTermSpecificity() == Ribonucleotide::FIVE_PRIME)
        {
          if (!seq.hasFivePrimeMod()) { seq.setFivePrimeMod(f); }
        }
        else if (f->getTermSpecificity() == Ribonucleotide::THREE_PRIME)
        {
          if (!seq.hasThreePrimeMod()) { seq.setThreePrimeMod(f); }
        }
      }
    );

    // iterate over each nucleotide
    size_t residue_index(0);
    for (auto const & r : seq)
    {
      // skip already modified residue
      if (r.isModified()) { ++residue_index; continue; }

      //set fixed modifications
      std::for_each(fixed_mods.begin(), fixed_mods.end(), [&seq, &residue_index, r] (ConstRibonucleotidePtr const & f)
        {
          // check if modification and current ribo match
          const std::string& code = r.getCode();
          if (code.size() == 1 && code[0] == f->getOrigin())
          {
            // replace the nucleoside with the modified version (skip 5'/3' modifications)
            if (f->getTermSpecificity() == Ribonucleotide::ANYWHERE)
            {
              seq.set(residue_index, f);
            }
          }
        }
      );
      ++residue_index;
    }
  }

  map<int, vector<ModifiedNASequenceGenerator::ConstRibonucleotidePtr>>
  ModifiedNASequenceGenerator::compatibleModifications_(const set<ConstRibonucleotidePtr>& var_mods, const NASequence& seq)
  {
    map<int, vector<ConstRibonucleotidePtr>> map_compatibility;

    const int FIVE_PRIME_MODIFICATION_INDEX = -1;
    const int THREE_PRIME_MODIFICATION_INDEX = -2;

    // set terminal modifications, if any are specified
    std::for_each(var_mods.begin(), var_mods.end(), [&seq, &map_compatibility, &FIVE_PRIME_MODIFICATION_INDEX, &THREE_PRIME_MODIFICATION_INDEX] (ConstRibonucleotidePtr const & v)
      {
        if (v->getTermSpecificity() == Ribonucleotide::FIVE_PRIME)
        {
          if (!seq.hasFivePrimeMod()) { map_compatibility[FIVE_PRIME_MODIFICATION_INDEX].push_back(v); }
        }
        else if (v->getTermSpecificity() == Ribonucleotide::THREE_PRIME)
        {
          if (!seq.hasThreePrimeMod()) { map_compatibility[THREE_PRIME_MODIFICATION_INDEX].push_back(v); }
        }
      });

    size_t residue_index(0);
    for (auto const & r : seq)
    {
      // skip already modified residues
      if (r.isModified())
      { 
        ++residue_index;
        continue;
      }

      //determine compatibility of variable modifications
      std::for_each(var_mods.begin(), var_mods.end(), [&residue_index, &r, &map_compatibility](ConstRibonucleotidePtr const & v)
      {
        // check if modification and current ribo match
        const std::string& code = r.getCode();
        if (code.size() == 1 && code[0] == v->getOrigin())
        {
          if (v->getTermSpecificity() == Ribonucleotide::ANYWHERE)
          {
            map_compatibility[static_cast<int>(residue_index)].push_back(v);
          }
        }
      });
      ++residue_index;
    }

    return map_compatibility;
  }

  // static
  void ModifiedNASequenceGenerator::applyVariableModifications(const set<ConstRibonucleotidePtr>& var_mods,
                                                               const NASequence& seq,
                                                               size_t max_variable_mods_per_seq,
                                                               vector<NASequence>& all_modified_seqs,
                                                               bool keep_unmodified)
  {
    // no variable modifications specified or no variable mods allowed? no compatibility map needs to be build
    if (var_mods.empty() || max_variable_mods_per_seq == 0)
    {
      // if unmodified seqs should be kept return the original list of digested seqs
      if (keep_unmodified) { all_modified_seqs.push_back(seq); }
      return;
    }

    // if there is at most one variable modification allowed for a seq we don't need combinatoric placement and can resort to a faster implementation
    if (max_variable_mods_per_seq == 1)
    {
      applyAtMostOneVariableModification_(var_mods, seq, all_modified_seqs, keep_unmodified);
      return;
    }

    // keep a list of all possible modifications of this seq
    vector<NASequence> modified_seqs;

    // only add unmodified version if flag is set (default)
    if (keep_unmodified) { modified_seqs.push_back(seq); }

    auto map_compatibility = compatibleModifications_(var_mods, seq);

    // Check if no compatible site that can be modified by variable
    // modification. If so just return seqs without variable modifications.
    const size_t & compatible_mod_sites = map_compatibility.size();
    if (compatible_mod_sites == 0)
    {
      if (keep_unmodified) { all_modified_seqs.push_back(seq); }
      return;
    }

    // generate powerset of max_variable_mods_per_seq sized subset of all compatible modification sites
    size_t max_placements = std::min(max_variable_mods_per_seq, compatible_mod_sites);
    for (size_t n_var_mods = 1; n_var_mods <= max_placements; ++n_var_mods)
    {
      // enumerate all modified seqs with n_var_mods variable modified residues
      size_t zeros = std::max((size_t)0, compatible_mod_sites - n_var_mods);
      vector<bool> subset_mask;

      for (size_t i = 0; i != compatible_mod_sites; ++i)
      {
        // create mask 000011 to select last (e.g. n_var_mods = 2) two compatible sites as subset from the set of all compatible sites
        if (i < zeros)
        {
          subset_mask.push_back(false);
        }
        else
        {
          subset_mask.push_back(true);
        }
      }

      // generate all subsets of compatible sites {000011, ... , 101000, 110000} with current number of allowed variable modifications per seq
      do
      {
        // create subset indices e.g.{4,12} from subset mask e.g. 1010000 corresponding to the positions in the seq sequence
        vector<int> subset_indices;
        map<int, vector<ConstRibonucleotidePtr>>::const_iterator mit = map_compatibility.begin();
        for (size_t i = 0; i != compatible_mod_sites; ++i, ++mit)
        {
          if (subset_mask[i])
          { 
            subset_indices.push_back(mit->first);
          }
        }

        // now enumerate all modifications
        recurseAndGenerateVariableModifiedSequences_(subset_indices, map_compatibility, 0, seq, modified_seqs);
      }
      while (next_permutation(subset_mask.begin(), subset_mask.end()));
    }
    // add modified version of the current seq to the list of all seqs
    all_modified_seqs.insert(all_modified_seqs.end(), modified_seqs.begin(), modified_seqs.end());
  }


  ModifiedNASequenceGenerator::PlacementStatistics
  ModifiedNASequenceGenerator::applyVariableModificationsConstrained(const set<ConstRibonucleotidePtr>& var_mods,
                                                                     const NASequence& seq,
                                                                     Size max_variable_mods,
                                                                     vector<NASequence>& output,
                                                                     const PlacementConstraints& constraints,
                                                                     bool keep_original)
  {
    if (! constraints.accepts || constraints.max_nodes == 0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Constrained placement needs a predicate and a positive node bound");
    }
    PlacementStatistics stats;
    const Size initial_size = output.size();
    const auto compatibility = max_variable_mods > 0 ? compatibleModifications_(var_mods, seq) : map<int, vector<ConstRibonucleotidePtr>> {};
    max_variable_mods = std::min(max_variable_mods, compatibility.size());
    const auto massOf
      = [&constraints](const NASequence& sequence) { return constraints.use_average_mass ? sequence.getAverageWeight() : sequence.getMonoWeight(); };
    const auto place = [](NASequence& sequence, int site, ConstRibonucleotidePtr modification) {
      if (site == -2) { sequence.setThreePrimeMod(modification); }
      else if (site == -1) { sequence.setFivePrimeMod(modification); }
      else
      {
        sequence.set(site, modification);
      }
    };
    const double base_mass = massOf(seq);
    map<int, vector<double>> deltas;
    map<int, std::pair<double, double>> ranges;
    bool finite_bounds = std::isfinite(base_mass);
    for (const auto& [site, choices] : compatibility)
    {
      auto& values = deltas[site];
      for (const auto* modification : choices)
      {
        // Exact sequence chemistry also handles terminal groups and linkage changes.
        NASequence modified = seq;
        place(modified, site, modification);
        const double delta = massOf(modified) - base_mass;
        finite_bounds = finite_bounds && std::isfinite(delta);
        values.push_back(delta);
      }
      const auto [low, high] = std::minmax_element(values.begin(), values.end());
      ranges.emplace(site, std::pair(*low, *high));
    }
    const auto alternativesOf = [](const NASequence& sequence) {
      ResidueAlternatives result;
      result.reserve(sequence.size());
      for (Size i = 0; i < sequence.size(); ++i)
      {
        result.push_back({sequence[i]});
      }
      return result;
    };
    const auto accepts = [&](double low, double high, const ResidueAlternatives& alternatives) {
      if (++stats.visited_nodes > constraints.max_nodes)
      {
        stats.enumeration_fallback = true;
        return false;
      }
      if (! finite_bounds)
      {
        low = -std::numeric_limits<double>::infinity();
        high = std::numeric_limits<double>::infinity();
      }
      else
      {
        // Outward slack covers floating-point regrouping of formula differences.
        const double slack = 1e-8 * std::max<Size>(1, seq.size() + compatibility.size())
                             + 64 * std::numeric_limits<double>::epsilon() * std::max(std::abs(low), std::abs(high));
        low -= slack;
        high += slack;
      }
      if (constraints.accepts(low, high, alternatives)) { return true; }
      ++stats.pruned_branches;
      return false;
    };

    // A budget-aware outer bound can reject the entire oligo before selecting sites.
    auto root_alternatives = alternativesOf(seq);
    vector<double> lower, upper;
    if (max_variable_mods > 0)
    {
      for (const auto& [site, choices] : compatibility)
      {
        lower.push_back(std::min(0.0, ranges.at(site).first));
        upper.push_back(std::max(0.0, ranges.at(site).second));
        if (site >= 0) { root_alternatives[site].insert(root_alternatives[site].end(), choices.begin(), choices.end()); }
      }
    }
    std::sort(lower.begin(), lower.end());
    std::sort(upper.begin(), upper.end(), std::greater<double>());
    double low = base_mass, high = base_mass;
    for (Size i = 0; i < max_variable_mods; ++i)
    {
      low += lower[i];
      high += upper[i];
    }
    if (! accepts(low, high, root_alternatives)) { return stats; }
    if (max_variable_mods == 0)
    {
      if (keep_original) { output.push_back(seq); }
      stats.emitted_variants = output.size() - initial_size;
      return stats;
    }

    const auto visit = [&](auto&& self, const vector<int>& sites, Size depth, const NASequence& current, double mass) -> bool {
      auto alternatives = alternativesOf(current);
      double branch_low = mass, branch_high = mass;
      for (Size i = depth; i < sites.size(); ++i)
      {
        const int site = sites[i];
        branch_low += ranges.at(site).first;
        branch_high += ranges.at(site).second;
        if (site >= 0) { alternatives[site] = compatibility.at(site); }
      }
      if (! accepts(branch_low, branch_high, alternatives)) { return ! stats.enumeration_fallback; }
      if (depth == sites.size())
      {
        output.push_back(current);
        return true;
      }
      const int site = sites[depth];
      const auto& choices = compatibility.at(site);
      for (Size i = 0; i < choices.size(); ++i)
      {
        NASequence modified = current;
        place(modified, site, choices[i]);
        if (! self(self, sites, depth + 1, modified, mass + deltas.at(site)[i])) { return false; }
      }
      return true;
    };
    if (keep_original) { visit(visit, {}, 0, seq, base_mass); }
    for (Size count = 1; count <= max_variable_mods && ! stats.enumeration_fallback; ++count)
    {
      vector<bool> mask(compatibility.size(), false);
      std::fill(mask.end() - count, mask.end(), true);
      do
      {
        vector<int> sites;
        Size index = 0;
        for (const auto& [site, choices] : compatibility)
        {
          if (mask[index++]) { sites.push_back(site); }
        }
        if (! visit(visit, sites, 0, seq, base_mass)) { break; }
      } while (std::next_permutation(mask.begin(), mask.end()));
    }
    if (stats.enumeration_fallback)
    {
      output.resize(initial_size);
      applyVariableModifications(var_mods, seq, max_variable_mods, output, keep_original);
    }
    stats.emitted_variants = output.size() - initial_size;
    return stats;
  }

  // static
  void ModifiedNASequenceGenerator::recurseAndGenerateVariableModifiedSequences_(
    const vector<int>& subset_indices,
    const map<int, vector<ConstRibonucleotidePtr>>& map_compatibility,
    int depth,
    const NASequence& current_seq,
    vector<NASequence>& modified_seqs)
  {
    const int FIVE_PRIME_MODIFICATION_INDEX = -1;
    const int THREE_PRIME_MODIFICATION_INDEX = -2;
    // cout << depth << " " << subset_indices.size() << " " << current_seq.toString() << endl;

    // end of recursion. Add the modified seq and return
    if (depth == (int)subset_indices.size())
    {
      modified_seqs.push_back(current_seq);
      return;
    }

    // get modifications compatible to residue at current seq position
    const int current_index = subset_indices[depth];

    auto const pos_mod_it = map_compatibility.find(current_index);
    const vector<ConstRibonucleotidePtr>& mods = pos_mod_it->second; // we don't need to check for .end as entry is guaranteed to exist

    for (auto const & m : mods)
    {
      // copy seq and apply modification
      NASequence new_seq = current_seq;
      if (current_index == THREE_PRIME_MODIFICATION_INDEX)
      {
        new_seq.setThreePrimeMod(m);
      }
      else if (current_index == FIVE_PRIME_MODIFICATION_INDEX)
      {
        new_seq.setFivePrimeMod(m);
      }
      else
      {
        new_seq.set(current_index, m);
      }

      // recurse with modified seq
      recurseAndGenerateVariableModifiedSequences_(subset_indices, map_compatibility, depth + 1, new_seq, modified_seqs);
    }
  }

  // static
  void ModifiedNASequenceGenerator::applyAtMostOneVariableModification_(
    const set<ConstRibonucleotidePtr>& var_mods,
    const NASequence& seq,
    vector<NASequence>& all_modified_seqs,
    bool keep_unmodified)
  {
    if (keep_unmodified)
    { 
      all_modified_seqs.push_back(seq);
    }

    for (const auto* mod : var_mods)
    {
      if (mod->getTermSpecificity() == Ribonucleotide::FIVE_PRIME && !seq.hasFivePrimeMod())
      {
        NASequence modified = seq;
        modified.setFivePrimeMod(mod);
        all_modified_seqs.push_back(modified);
      }
      else if (mod->getTermSpecificity() == Ribonucleotide::THREE_PRIME && !seq.hasThreePrimeMod())
      {
        NASequence modified = seq;
        modified.setThreePrimeMod(mod);
        all_modified_seqs.push_back(modified);
      }
    }
    // Walk backwards without constructing an iterator before begin().
    for (Size remaining = seq.size(); remaining > 0; --remaining)
    {
      const auto ribo_it = seq.cbegin() + (remaining - 1);
      // skip already modified residues
      if (ribo_it->isModified())
      { 
        continue;
      }
      size_t residue_index = ribo_it - seq.cbegin();

      // matches every variable modification to every site and return the new sequence with single modification
      std::for_each(var_mods.begin(), var_mods.end(),
                    [ribo_it, residue_index, &all_modified_seqs, &seq](ConstRibonucleotidePtr const & v)
      {
        // check if modification and current ribo match
        const std::string& code = ribo_it->getCode();
        if (v->getTermSpecificity() == Ribonucleotide::ANYWHERE && code.size() == 1 && code[0] == v->getOrigin())
        {
          NASequence new_seq = seq;
          new_seq.set(residue_index, v);
          all_modified_seqs.push_back(new_seq);
        }
      });
    }
  }
}

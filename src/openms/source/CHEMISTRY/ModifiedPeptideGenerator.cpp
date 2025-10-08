// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <boost/math/special_functions/binomial.hpp>

using std::vector;
using std::pair;
using std::map;

namespace OpenMS
{

  constexpr int ModifiedPeptideGenerator::N_TERM_MODIFICATION_INDEX = -1; // magic constant to distinguish N_TERM only modifications from ANYWHERE modifications placed at N-term residue
  constexpr int ModifiedPeptideGenerator::C_TERM_MODIFICATION_INDEX = -2; // magic constant to distinguish C_TERM only modifications from ANYWHERE modifications placed at C-term residue

  // static
  ModifiedPeptideGenerator::MapToResidueType ModifiedPeptideGenerator::getModifications(const StringList& modNames)
  {   
    vector<const ResidueModification*> modifications;
    
    for (const String& modification : modNames)
    {
      const ResidueModification* rm = ModificationsDB::getInstance()->getModification(modification);
      modifications.push_back(rm);
    }
    std::sort(modifications.begin(), modifications.end());
    return createResidueModificationToResidueMap_(modifications);
  }


  // static
  ModifiedPeptideGenerator::MapToResidueType ModifiedPeptideGenerator::createResidueModificationToResidueMap_(const vector<const ResidueModification*>& mods)
  {
    // create a lookup structure from ResidueModification (e.g., "Oxidation (M)" to the modified Residue* in ResidueDB"
    ModifiedPeptideGenerator::MapToResidueType m;
    for (auto const & r : mods)
    {
      String name = r->getFullId();
      bool is_terminal = r->getTermSpecificity() == ResidueModification::N_TERM || r->getTermSpecificity() == ResidueModification::C_TERM || r->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM || r->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM;
      if (!is_terminal)
      {
        auto residue = ResidueDB::getInstance()->getResidue(r->getOrigin());
        m.val[r] = ResidueDB::getInstance()->getModifiedResidue(residue, name);
      }
      else  // terminal modification
      {
        if (r->getOrigin() == 'X')
        { // no residue associated with strictly terminal modification
          m.val[r] = nullptr; 
        }
        else
        { // specific residue associated with strictly terminal modification
          auto residue = ResidueDB::getInstance()->getResidue(r->getOrigin());
          m.val[r] = ResidueDB::getInstance()->getModifiedResidue(residue, name);          
        }                
      }      
    }
    return m;
  }

  // static
  void ModifiedPeptideGenerator::applyFixedModifications(
    const MapToResidueType& fixed_mods, 
    AASequence& peptide)
  {
    // set terminal modifications for modifications without amino acid preference
    for (auto const& mr : fixed_mods.val)
    {
      const ResidueModification* f = mr.first;
      if (f->getTermSpecificity() == ResidueModification::N_TERM)
      {
        if (!peptide.hasNTerminalModification())
        {
          peptide.setNTerminalModification(f);
        }
      }
      else if (f->getTermSpecificity() == ResidueModification::C_TERM)
      {
        if (!peptide.hasCTerminalModification())
        {
          peptide.setCTerminalModification(f);
        }
      }
    }

    //iterate over each residue
    for (auto residue_it = peptide.begin(); residue_it != peptide.end(); ++residue_it)
    {
      // skip already modified residue
      if (residue_it->isModified())
      {
        continue;
      }
      Size residue_index = residue_it - peptide.begin();
      //set fixed modifications
      for (auto const& mr : fixed_mods.val)
      {
        const ResidueModification* f = mr.first;

        // check if amino acid match between modification and current residue
        if (residue_it->getOneLetterCode()[0] != f->getOrigin()) { continue; }

        // Term specificity is ANYWHERE on the peptide, C_TERM or N_TERM (currently no explicit support in OpenMS for protein C-term and protein N-term)
        const ResidueModification::TermSpecificity& term_spec = f->getTermSpecificity();
        if (term_spec == ResidueModification::ANYWHERE)
        {
          const Residue* r = mr.second; // map modification to the modified residue
          peptide.setModification(residue_index, r);
        }
        else if (term_spec == ResidueModification::C_TERM && residue_index == (peptide.size() - 1))
        {
          peptide.setCTerminalModification(f);
        }
        else if (term_spec == ResidueModification::N_TERM && residue_index == 0)
        {
          peptide.setNTerminalModification(f);
        }
      }
    }
  }

  // static
  void ModifiedPeptideGenerator::applyVariableModifications(
    const MapToResidueType& var_mods,
    const AASequence& peptide, 
    Size max_variable_mods_per_peptide, 
    vector<AASequence>& all_modified_peptides, 
    bool keep_unmodified)
  {
    // no variable modifications specified or no variable mods allowed? no compatibility map needs to be build
    if (var_mods.val.empty() || max_variable_mods_per_peptide == 0)
    {
      // if unmodified peptides should be kept return the original list of digested peptides
      if (keep_unmodified) { all_modified_peptides.push_back(peptide); }
      return;
    }

    // if there is at most one variable modification allowed for a peptide we don't need combinatoric placement and can reside to a faster implementation
    if (max_variable_mods_per_peptide == 1)
    {
      applyAtMostOneVariableModification_(
        var_mods,
        peptide,
        all_modified_peptides,
        keep_unmodified);
      return;
    }

    // iterate over each residue and build compatibility mapping describing
    // which amino acid (peptide index) is compatible with which modification
    map<int, vector<const ResidueModification*> > mod_compatibility;

    // set terminal modifications for modifications without amino acid preference
    for (auto const& mr : var_mods.val)
    {
      const ResidueModification* v = mr.first;

      if (v->getTermSpecificity() == ResidueModification::N_TERM)
      {
        if (!peptide.hasNTerminalModification())
        {
          mod_compatibility[N_TERM_MODIFICATION_INDEX].push_back(v);
        }
      }
      else if (v->getTermSpecificity() == ResidueModification::C_TERM)
      {
        if (!peptide.hasCTerminalModification())
        {
          mod_compatibility[C_TERM_MODIFICATION_INDEX].push_back(v);
        }
      }
    }

    
    for (auto residue_it = peptide.begin(); residue_it != peptide.end(); ++residue_it)
    {
      // skip already modified residues
      if (residue_it->isModified())
      {
        continue;
      }

      Size residue_index = residue_it - peptide.begin();

      //determine compatibility of variable modifications
      for (auto const& mr : var_mods.val)
      {
        const ResidueModification* v = mr.first;

        // check if amino acid match between modification and current residue
        if (residue_it->getOneLetterCode()[0] != v->getOrigin())
        {
          continue;
        }

        // Term specificity is ANYWHERE on the peptide, C_TERM or N_TERM
        // (currently no explicit support in OpenMS for protein C-term and
        // protein N-term)
        // TODO This is not true anymore!
        const ResidueModification::TermSpecificity& term_spec = v->getTermSpecificity();
        if (term_spec == ResidueModification::ANYWHERE)
        {
          mod_compatibility[residue_index].push_back(v);
        }
        // TODO think about if it really is the same case as the one above.
        else if (term_spec == ResidueModification::C_TERM && residue_index == (peptide.size() - 1))
        {
          mod_compatibility[C_TERM_MODIFICATION_INDEX].push_back(v);
        }
        else if (term_spec == ResidueModification::N_TERM && residue_index == 0)
        {
          mod_compatibility[N_TERM_MODIFICATION_INDEX].push_back(v);
        }
      }
    }

    Size max_placements = std::min(max_variable_mods_per_peptide, mod_compatibility.size());

    // stores all variants with how many modifications they already have
    vector<pair<size_t, vector<AASequence>>> mod_peps_w_depth = {{0, {peptide}}};
    Size num_res = 0;
    for (Size s(0); s <= max_placements; ++s)
    {
      num_res += boost::math::binomial_coefficient<double>(mod_compatibility.size(), s);
    }
    mod_peps_w_depth.reserve(num_res);
    auto rit = mod_compatibility.rbegin();
    for (; rit != mod_compatibility.rend(); ++rit)
    {
      const auto& idx = rit->first;
      const auto& mods = rit->second;
      // copy the complete sequences from last iteration
      auto tmp = mod_peps_w_depth;
      for (auto& [old_depth, old_variants] : tmp)
      {
        // extends mod_peps_w_depth by adding variants with the next mod, if max_placements is not reached
        if (old_depth < max_placements)
        {
          applyAllModsAtIdxAndExtend_(old_variants, idx, mods, var_mods);
          mod_peps_w_depth.emplace_back(old_depth + 1, std::move(old_variants));
        }
      }
    }
    
    // move sequences from mod_peps_w_depth into result. Skip the initial peptide if desired.
    for (auto& [depth, seqs] : mod_peps_w_depth)
    {
      if (depth != 0 || keep_unmodified)
      {
        all_modified_peptides.insert(
          all_modified_peptides.end(), 
          make_move_iterator(seqs.begin()), 
          make_move_iterator(seqs.end())); 
      }
    }
  }

  // static
  void ModifiedPeptideGenerator::applyAtMostOneVariableModification_(
    const MapToResidueType& var_mods, 
    const AASequence& peptide, 
    vector<AASequence>& all_modified_peptides, 
    bool keep_unmodified)
  {
    if (keep_unmodified)
    {
      all_modified_peptides.push_back(peptide);
    }

    // we want the same behavior as for the slower function... we would need a reverse iterator here that AASequence doesn't provide
    for (auto residue_it = peptide.end() - 1; residue_it != peptide.begin() - 1; --residue_it)
    {
      // skip already modified residues
      if (residue_it->isModified())
      {
        continue;
      }

      Size residue_index = residue_it - peptide.begin();

      // determine compatibility of variable modifications
      for (auto const & mr : var_mods.val)
      {
        const ResidueModification* v = mr.first;

        const char r = residue_it->getOneLetterCode()[0];
        // check if amino acid match between modification and current residue
        if (r != v->getOrigin()) { continue; }

        // Term specificity is ANYWHERE on the peptide, C_TERM or N_TERM (currently no explicit support in OpenMS for protein C-term and protein N-term)
        const ResidueModification::TermSpecificity& term_spec = v->getTermSpecificity();
        bool is_compatible(false);
        if (term_spec == ResidueModification::ANYWHERE)
        {
          is_compatible = true;
        }
        else if (term_spec == ResidueModification::C_TERM && residue_index == (peptide.size() - 1))
        {
          is_compatible = true;
        }
        else if (term_spec == ResidueModification::N_TERM && residue_index == 0)
        {
          is_compatible = true;
        }

        // residue modification an be placed at current position? Then generate modified peptide.
        if (is_compatible)
        {
          AASequence new_peptide = peptide;
          new_peptide.setModification(residue_index, mr.second); // set modified Residue          
          all_modified_peptides.push_back(std::move(new_peptide));
        }
      }
    }
  }


  void ModifiedPeptideGenerator::applyAllModsAtIdxAndExtend_(vector<AASequence>& original_sequences, int idx_to_modify, const vector<const ResidueModification*>& mods, const ModifiedPeptideGenerator::MapToResidueType& var_mods)
  {
    Size end = original_sequences.size();
    original_sequences.reserve(end * mods.size());
    for (Size s(1); s < mods.size(); ++s)
    {
      original_sequences.insert(original_sequences.end(),original_sequences.begin(), original_sequences.begin()+end);
    }
    for (Size cnt(0); cnt < mods.size(); ++cnt) // apply first mod later
    {
      for (Size i(0); i < end; i++)
      {
        applyModToPep_(original_sequences[cnt * end + i], idx_to_modify, mods[cnt], var_mods);
      }
    }
  }


  void ModifiedPeptideGenerator::applyModToPep_(AASequence& current_peptide, int current_index, const ResidueModification* m, const ModifiedPeptideGenerator::MapToResidueType& var_mods)
  {
      if (current_index == C_TERM_MODIFICATION_INDEX)
      {
        current_peptide.setCTerminalModification(m);
      }
      else if (current_index == N_TERM_MODIFICATION_INDEX)
      {
        current_peptide.setNTerminalModification(m);
      }
      else
      {
        const Residue* r = var_mods.val.at(m); // map modification to the modified residue
        current_peptide.setModification(current_index, r); // set modified Residue       
      }
  }
}


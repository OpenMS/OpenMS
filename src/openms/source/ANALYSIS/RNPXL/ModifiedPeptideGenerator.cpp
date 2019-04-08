// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
using std::vector;
using std::map;

namespace OpenMS
{

 // static
 ModifiedPeptideGenerator::MapToResidueType ModifiedPeptideGenerator::getModifications(const StringList& modNames)
 {   
   vector<const ResidueModification*> modifications;
  
   for (String modification : modNames)
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
    // create a lookup structure from ResidueModification (e.g., "Oxidiation (M)" to the modified Residue* in ResidueDB"
    ModifiedPeptideGenerator::MapToResidueType m;
    for (auto const & r : mods)
    {
      String name = r->getFullId();
      bool is_terminal = r->getTermSpecificity() == ResidueModification::N_TERM || r->getTermSpecificity() == ResidueModification::C_TERM;

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

    const int N_TERM_MODIFICATION_INDEX = -1; // magic constant to distinguish N_TERM only modifications from ANYWHERE modifications placed at N-term residue
    const int C_TERM_MODIFICATION_INDEX = -2; // magic constant to distinguish C_TERM only modifications from ANYWHERE modifications placed at C-term residue

    //keep a list of all possible modifications of this peptide
    vector<AASequence> modified_peptides;

    // only add unmodified version if flag is set (default)
    if (keep_unmodified)
    {
      modified_peptides.push_back(peptide);
    }

    // iterate over each residue and build compatibility mapping describing
    // which amino acid (peptide index) is compatible with which modification
    map<int, vector<const ResidueModification*> > map_compatibility;

    // set terminal modifications for modifications without amino acid preference
    for (auto const& mr : var_mods.val)
    {
      const ResidueModification* v = mr.first;

      if (v->getTermSpecificity() == ResidueModification::N_TERM)
      {
        if (!peptide.hasNTerminalModification())
        {
          map_compatibility[N_TERM_MODIFICATION_INDEX].push_back(v);
        }
      }
      else if (v->getTermSpecificity() == ResidueModification::C_TERM)
      {
        if (!peptide.hasCTerminalModification())
        {
          map_compatibility[C_TERM_MODIFICATION_INDEX].push_back(v);
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
        const ResidueModification::TermSpecificity& term_spec = v->getTermSpecificity();
        if (term_spec == ResidueModification::ANYWHERE)
        {
          map_compatibility[static_cast<int>(residue_index)].push_back(v);
        }
        else if (term_spec == ResidueModification::C_TERM && residue_index == (peptide.size() - 1))
        {
          map_compatibility[C_TERM_MODIFICATION_INDEX].push_back(v);
        }
        else if (term_spec == ResidueModification::N_TERM && residue_index == 0)
        {
          map_compatibility[N_TERM_MODIFICATION_INDEX].push_back(v);
        }
      }
    }

    // Check if no compatible site that can be modified by variable
    // modification. If so just return peptides without variable modifications.
    const Size compatible_mod_sites = map_compatibility.size();
    if (compatible_mod_sites == 0)
    {
      if (keep_unmodified)
      {
        all_modified_peptides.push_back(peptide);
      }
      return;
    }

    // generate powerset of max_variable_mods_per_peptide sized subset of all compatible modification sites
    Size max_placements = std::min(max_variable_mods_per_peptide, compatible_mod_sites);
    for (Size n_var_mods = 1; n_var_mods <= max_placements; ++n_var_mods)
    {
      // enumerate all modified peptides with n_var_mods variable modified residues
      Size zeros = std::max((Size)0, compatible_mod_sites - n_var_mods);
      vector<bool> subset_mask;

      for (Size i = 0; i != compatible_mod_sites; ++i)
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

      // generate all subsets of compatible sites {000011, ... , 101000, 110000} with current number of allowed variable modifications per peptide
      do
      {
        // create subset indices e.g.{4,12} from subset mask e.g. 1010000 corresponding to the positions in the peptide sequence
        vector<int> subset_indices;
        map<int, vector<const ResidueModification*> >::const_iterator mit = map_compatibility.begin();
        for (Size i = 0; i != compatible_mod_sites; ++i, ++mit)
        {
          if (subset_mask[i])
          {
            subset_indices.push_back(mit->first);
          }
        }

        // now enumerate all modifications
        recurseAndGenerateVariableModifiedPeptides_(subset_indices, map_compatibility, var_mods, 0, peptide, modified_peptides);
      } while (next_permutation(subset_mask.begin(), subset_mask.end()));
    }
    // add modified version of the current peptide to the list of all peptides
    
    all_modified_peptides.insert(
      all_modified_peptides.end(), 
      make_move_iterator(modified_peptides.begin()), 
      make_move_iterator(modified_peptides.end())); 
      
  }


  // static
  void ModifiedPeptideGenerator::recurseAndGenerateVariableModifiedPeptides_(
    const vector<int>& subset_indices, 
    const map<int, vector<const ResidueModification*> >& map_compatibility, 
    const MapToResidueType& var_mods, 
    int depth, 
    const AASequence& current_peptide, 
    vector<AASequence>& modified_peptides)
  {
    const int N_TERM_MODIFICATION_INDEX = -1; // magic constant to distinguish N_TERM only modifications from ANYWHERE modifications placed at N-term residue
    const int C_TERM_MODIFICATION_INDEX = -2; // magic constant to distinguish C_TERM only modifications from ANYWHERE modifications placed at C-term residue

    // cout << depth << " " << subset_indices.size() << " " << current_peptide.toString() << endl;

    // end of recursion. Add the modified peptide and return
    if (depth == (int)subset_indices.size())
    {
      modified_peptides.push_back(current_peptide);
      return;
    }

    // get modifications compatible to residue at current peptide position
    const int current_index = subset_indices[depth];

    map<int, vector<const ResidueModification*> >::const_iterator pos_mod_it = map_compatibility.find(current_index);
    const vector<const ResidueModification*>& mods = pos_mod_it->second; // we don't need to check for .end as entry is guaranteed to exist

    for (const ResidueModification* m : mods)
    {

      // copy peptide and apply modification
      AASequence new_peptide = current_peptide;
      if (current_index == C_TERM_MODIFICATION_INDEX)
      {
        new_peptide.setCTerminalModification(m);
      }
      else if (current_index == N_TERM_MODIFICATION_INDEX)
      {
        new_peptide.setNTerminalModification(m);
      }
      else
      {
        const Residue* r = var_mods.val.at(m); // map modification to the modified residue
        new_peptide.setModification(current_index, r); // set modified Residue          
      }

      // recurse with modified peptide
      recurseAndGenerateVariableModifiedPeptides_(subset_indices, map_compatibility, var_mods, depth + 1, new_peptide, modified_peptides);
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
}


// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

using namespace std;

#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
namespace OpenMS
{
  // static
  void ModifiedPeptideGenerator::applyFixedModifications(const vector<ResidueModification>::const_iterator& fixed_mods_begin, const vector<ResidueModification>::const_iterator& fixed_mods_end, vector<AASequence>::iterator digested_peptides_begin, vector<AASequence>::iterator digested_peptides_end)
  {
    for (vector<AASequence>::iterator peptide_it = digested_peptides_begin; peptide_it != digested_peptides_end; ++peptide_it)
    {
      AASequence& peptide_seq = *peptide_it;
      //iterate over each residue
      for (AASequence::ConstIterator residue_it = peptide_seq.begin(); residue_it != peptide_seq.end(); ++residue_it)
      {
        // skip already modified residue
        if (residue_it->isModified())
        {
          continue;
        }
        Size residue_index = residue_it - peptide_seq.begin();
        //set fixed modifications
        for (vector<ResidueModification>::const_iterator fixed_it = fixed_mods_begin; fixed_it != fixed_mods_end; ++fixed_it)
        {
          // check if amino acid match between modification and current residue
          if (residue_it->getOneLetterCode() != fixed_it->getOrigin())
          {
            continue;
          }

          // Term specifity is ANYWHERE on the peptide, C_TERM or N_TERM (currently no explicit support in OpenMS for protein C-term and protein N-term)
          const ResidueModification::Term_Specificity& term_spec = fixed_it->getTermSpecificity();
          if (term_spec == ResidueModification::ANYWHERE)
          {
            peptide_seq.setModification(residue_index, fixed_it->getFullName());
          } else if (term_spec == ResidueModification::C_TERM && residue_index == (peptide_seq.size() - 1))
          {
            peptide_seq.setCTerminalModification(fixed_it->getFullName());
          } else if (term_spec == ResidueModification::N_TERM && residue_index == 0)
          {
            peptide_seq.setNTerminalModification(fixed_it->getFullName());
          }
        }
      }
    }
  }

  // static
  void ModifiedPeptideGenerator::applyVariableModifications(const vector<ResidueModification>::const_iterator& var_mods_begin, const vector<ResidueModification>::const_iterator& var_mods_end, vector<AASequence>::iterator digested_peptides_begin, vector<AASequence>::iterator digested_peptides_end, Size max_variable_mods_per_peptide, vector<AASequence>& all_modified_peptides, bool keep_unmodified)
  {
    // no variable modifications specified? no compatibility map needs to be build
    if (var_mods_begin == var_mods_end)
    {
      // if unmodifed peptides should be kept return the original list of digested peptides
      if (keep_unmodified)
      {
        all_modified_peptides.insert(all_modified_peptides.end(), digested_peptides_begin, digested_peptides_end);
      }
      return;
    }

    const int N_TERM_MODIFICATION_INDEX = -1; // magic constant to distinguish N_TERM only modifications from ANYWHERE modifications placed at N-term residue
    const int C_TERM_MODIFICATION_INDEX = -2; // magic constant to distinguish C_TERM only modifications from ANYWHERE modifications placed at C-term residue

    for (vector<AASequence>::const_iterator peptide_it = digested_peptides_begin; peptide_it != digested_peptides_end; ++peptide_it)
    {
      const AASequence& peptide_seq = *peptide_it;

      //keep a list of all possible modifications of this peptide
      vector<AASequence> modified_peptides;

      // only add unmodified version if flag is set (default)
      if (keep_unmodified)
      {
        modified_peptides.push_back(peptide_seq);
      }

      //iterate over each residue and build compatibility mapping describing which amino acid (peptide index) is compatible with which modification
      map<int, vector<ResidueModification> > map_compatibility;

      for (AASequence::ConstIterator residue_it = peptide_seq.begin(); residue_it != peptide_seq.end(); ++residue_it)
      {
        // skip already modified residues
        if (residue_it->isModified())
        {
          continue;
        }

        Size residue_index = residue_it - peptide_seq.begin();

        //determine compatibility of variable modifications
        for (vector<ResidueModification>::const_iterator variable_it = var_mods_begin; variable_it != var_mods_end; ++variable_it)
        {
          // check if amino acid match between modification and current residue
          if (residue_it->getOneLetterCode() != variable_it->getOrigin())
          {
            continue;
          }

          // Term specifity is ANYWHERE on the peptide, C_TERM or N_TERM (currently no explicit support in OpenMS for protein C-term and protein N-term)
          const ResidueModification::Term_Specificity& term_spec = variable_it->getTermSpecificity();
          if (term_spec == ResidueModification::ANYWHERE)
          {
            map_compatibility[residue_index].push_back(*variable_it);
          } else if (term_spec == ResidueModification::C_TERM && residue_index == (peptide_seq.size() - 1))
          {
            map_compatibility[C_TERM_MODIFICATION_INDEX].push_back(*variable_it);
          } else if (term_spec == ResidueModification::N_TERM && residue_index == 0)
          {
            map_compatibility[N_TERM_MODIFICATION_INDEX].push_back(*variable_it);
          }
        }
      }

      // Check if no compatible site that can be modified by variable modification. If so just return peptides without variable modifications.
      const Size compatible_mod_sites = map_compatibility.size();
      if (compatible_mod_sites == 0)
      {
        if (keep_unmodified)
        {
          all_modified_peptides.push_back(peptide_seq);
        }
        continue;
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
          } else
          {
            subset_mask.push_back(true);
          }
        }

        // generate all subsets of compatible sites {000011, ... , 101000, 110000} with current number of allowed variable modifications per peptide
        do
        {
          // create subset indices e.g.{4,12} from subset mask e.g. 1010000 corresponding to the positions in the peptide sequence
          vector<int> subset_indices;
          map<int, vector<ResidueModification> >::const_iterator mit = map_compatibility.begin();
          for (Size i = 0; i != compatible_mod_sites; ++i, ++mit)
          {
            if (subset_mask[i])
            {
              subset_indices.push_back(mit->first);
            }
          }

          // now enumerate all modifications
          recurseAndGenerateVariableModifiedPeptides_(subset_indices, map_compatibility, 0, peptide_seq, modified_peptides);
        } while (next_permutation(subset_mask.begin(), subset_mask.end()));
      }
      // add modified version of the current peptide to the list of all peptides
      all_modified_peptides.insert(all_modified_peptides.end(), modified_peptides.begin(), modified_peptides.end());
    } // for current peptide
  }

  // static
  void ModifiedPeptideGenerator::recurseAndGenerateVariableModifiedPeptides_(const vector<int>& subset_indices, const map<int, vector<ResidueModification> >& map_compatibility, int depth, AASequence current_peptide, vector<AASequence>& modified_peptides)
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
    const vector<ResidueModification>& mods = map_compatibility.at(current_index);

    for (vector<ResidueModification>::const_iterator mod_it = mods.begin(); mod_it != mods.end(); ++mod_it)
    {
      // copy peptide and apply modification
      AASequence new_peptide = current_peptide;
      if (current_index == C_TERM_MODIFICATION_INDEX)
      {
        new_peptide.setCTerminalModification(mod_it->getFullName());
      }
      else if (current_index == N_TERM_MODIFICATION_INDEX)
      {
        new_peptide.setNTerminalModification(mod_it->getFullName());
      } else
      {
        new_peptide.setModification(current_index, mod_it->getFullName());
      }

      // recurse with modified peptide
      recurseAndGenerateVariableModifiedPeptides_(subset_indices, map_compatibility, depth + 1, new_peptide, modified_peptides);
    }
  }

  //static
}


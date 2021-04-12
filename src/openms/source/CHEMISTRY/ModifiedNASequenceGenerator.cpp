// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/CHEMISTRY/ModifiedNASequenceGenerator.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CHEMISTRY/NASequence.h>

#include <unordered_map>

using namespace std;

namespace OpenMS
{
  // static
  void ModifiedNASequenceGenerator::applyFixedModifications(
    const set<ConstRibonucleotidePtr>& fixed_mods,
    NASequence& seq)
  {
    // apply modifications at chain ends
    std::for_each(fixed_mods.begin(), fixed_mods.end(), [&seq](ConstRibonucleotidePtr f)
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
          const String& code = r.getCode();
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


  void ModifiedNASequenceGenerator::addModToSequences_(
    vector<ModSeqInfo>& temp_seqs, Size n_temp_seqs,
    vector<NASequence>& finished_seqs,
    const function<bool(NASequence&, Int&)>& applyMod)
  {
    for (Size i = 0; i < n_temp_seqs; ++i)
    {
      NASequence new_seq = temp_seqs[i].seq;
      Int missed_cleavages_left = temp_seqs[i].missed_cleavages_left;
      bool success = applyMod(new_seq, missed_cleavages_left);
      if (!success) continue; // no missed cleavages left, can't add inosine
      Size var_mods_left = temp_seqs[i].var_mods_left - 1;
      if (var_mods_left > 0)
      {
        temp_seqs.emplace_back(new_seq, var_mods_left, missed_cleavages_left);
      }
      else
      {
        finished_seqs.push_back(new_seq);
      }
    }
  }


  void ModifiedNASequenceGenerator::applyVariableModifications(
      const set<ConstRibonucleotidePtr>& var_mods, const NASequence& seq,
      Size max_var_mods, vector<NASequence>& all_modified_seqs,
      bool keep_unmodified, Int max_missed_cleavages)
  {
    all_modified_seqs.clear();
    if (keep_unmodified) all_modified_seqs.push_back(seq);
    if (var_mods.empty() || (max_var_mods == 0))
    {
      return;
    }

    // generate residue/mod. compatibility map:
    unordered_map<Size, set<ConstRibonucleotidePtr>> compatible_mods;
    set<ConstRibonucleotidePtr> compatible_5p_mods, compatible_3p_mods;
    for (ConstRibonucleotidePtr mod : var_mods)
    {
      switch (mod->getTermSpecificity())
      {
      case Ribonucleotide::FIVE_PRIME:
        if (!seq.hasFivePrimeMod()) compatible_5p_mods.insert(mod);
        break;
      case Ribonucleotide::THREE_PRIME:
        if (!seq.hasThreePrimeMod()) compatible_3p_mods.insert(mod);
        break;
      default: // case Ribonucleotide::ANYWHERE
        char origin = mod->getOrigin();
        for (Size i = 0; i < seq.size(); ++i)
        {
          const String& code = seq[i]->getCode();
          if ((code.size() == 1) && (code[0] == origin))
          {
            compatible_mods[i].insert(mod);
          }
        }
      }
    }

    // stop if there aren't any possible mod. placements:
    if (compatible_mods.empty() && compatible_5p_mods.empty() && (compatible_3p_mods.empty()))
    {
      return;
    }

    // buffer of sequences that can accept further mods. (and how many):
    vector<ModSeqInfo> temp_seqs;
    // starting with the original sequence, add one (more) possible mod. each time:
    temp_seqs.emplace_back(seq, max_var_mods, max_missed_cleavages);
    for (ConstRibonucleotidePtr mod : compatible_5p_mods)
    {
      addModToSequences_(temp_seqs, 1, all_modified_seqs,
                         [mod](NASequence& new_seq, Int& /* ignored */)
                         {
                           new_seq.setFivePrimeMod(mod);
                           return true;
                         });
    }
    for (const auto& comp_pair : compatible_mods)
    {
      Size pos = comp_pair.first;
      // only apply mods to sequences that are already in the buffer now:
      Size n_temp_seqs = temp_seqs.size();
      for (ConstRibonucleotidePtr mod : comp_pair.second)
      {
        if ((max_missed_cleavages < 0) || (mod->getCode() != "I"))
        {
          addModToSequences_(temp_seqs, n_temp_seqs, all_modified_seqs,
                             [mod, pos](NASequence& new_seq, Int& /* ignored */)
                             {
                               new_seq[pos] = mod;
                               return true;
                             });
        }
        else // special case for inosines (add RNase T1 cleavage sites)
        {
          addModToSequences_(temp_seqs, n_temp_seqs, all_modified_seqs,
                             [mod, pos](NASequence& new_seq, Int& missed_cleavages_left)
                             {
                               if (missed_cleavages_left <= 0) return false;
                               new_seq[pos] = mod;
                               --missed_cleavages_left;
                               return true;
                             });

        }
      }
    }
    Size n_temp_seqs = temp_seqs.size();
    for (ConstRibonucleotidePtr mod : compatible_3p_mods)
    {
      addModToSequences_(temp_seqs, n_temp_seqs, all_modified_seqs,
                         [mod](NASequence& new_seq, Int& /* ignored */)
                         {
                           new_seq.setThreePrimeMod(mod);
                           return true;
                         });
    }

    // add "partially" modified sequences to the output:
    for (Size i = 1; i < temp_seqs.size(); ++i)
    {
      all_modified_seqs.push_back(temp_seqs[i].seq);
    }
  }
}

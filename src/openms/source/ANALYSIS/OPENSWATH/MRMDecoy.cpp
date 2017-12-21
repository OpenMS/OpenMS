// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <map>
#include <utility> //for pair
#include <string>
#include <vector>
#include <algorithm> // for sort

namespace OpenMS
{
  std::vector<std::pair<std::string::size_type, std::string> > MRMDecoy::find_all_tryptic(std::string sequence)
  {
    std::vector<std::pair<std::string::size_type, std::string> > idx;
    std::vector<std::string> pattern;
    pattern.push_back("K");
    pattern.push_back("R");
    pattern.push_back("P");

    for (Size i = 0; i < sequence.size(); i++)
    {
      for (Size j = 0; j < pattern.size(); j++)
      {
        if (sequence.substr(i, 1) == pattern[j])
        {
          std::pair<std::string::size_type, std::string> idx_pair(i, pattern[j]);
          idx.push_back(idx_pair);
        }
      }
    }
    return idx;
  }

  float MRMDecoy::AASequenceIdentity(const String& sequence, const String& decoy)
  {
    OPENMS_PRECONDITION(sequence.size() == decoy.size(), "Cannot compare two sequences of unequal length");

    std::vector<char> sequence_v(sequence.begin(), sequence.end());
    std::vector<char> decoy_v(decoy.begin(), decoy.end());
    int running = 0;
    for (Size i = 0; i < sequence_v.size(); i++)
    {
      if (sequence_v[i] == decoy_v[i])
      {
        running += 1;
      }
    }
    double identity = (double) running / sequence_v.size();
    return identity;
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::shufflePeptide(
    OpenMS::TargetedExperiment::Peptide peptide, double identity_threshold, int seed,
    int max_attempts, bool replace_aa_instead_append)
  {
#ifdef DEBUG_MRMDECOY
    std::cout << " shuffle peptide " << peptide.sequence << std::endl;
    seed = 41;
#endif
    if (seed == -1)
    {
      seed = time(nullptr);
    }
    OpenMS::TargetedExperiment::Peptide shuffled = peptide;

    boost::mt19937 generator(seed);
    boost::uniform_int<> uni_dist;
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > pseudoRNG(generator, uni_dist);

    typedef std::vector<std::pair<std::string::size_type, std::string> > IndexType;
    IndexType idx = MRMDecoy::find_all_tryptic(peptide.sequence);
    std::string aa[] =
    {
      "A", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "M", "F", "S", "T", "W",
      "Y", "V"
    };
    int aa_size = 17;

    int attempts = 0;
    // loop: copy the original peptide, attempt to shuffle it and check whether difference is large enough
    while (MRMDecoy::AASequenceIdentity(peptide.sequence, shuffled.sequence) > identity_threshold &&
           attempts < max_attempts)
    {
      shuffled = peptide;
      std::vector<Size> peptide_index;
      for (Size i = 0; i < peptide.sequence.size(); i++)
      {
        peptide_index.push_back(i);
      }

      // we erase the indices where K/P/R are (from the back / in reverse order
      // to not delete indices we access later)
      for (IndexType::reverse_iterator it = idx.rbegin(); it != idx.rend(); ++it)
      {
        peptide_index.erase(peptide_index.begin() + it->first);
      }

      // shuffle the peptide index (without the K/P/R which we leave in place)
      // one could also use std::random_shuffle here but then the code becomes
      // untestable since the implementation of std::random_shuffle differs
      // between libc++ (llvm/mac-osx) and libstdc++ (gcc) and VS
      // see also https://code.google.com/p/chromium/issues/detail?id=358564
      // the actual code here for the shuffling is based on the implementation of
      // std::random_shuffle in libstdc++
      if (peptide_index.begin() != peptide_index.end())
      {
        for (std::vector<Size>::iterator pI_it = peptide_index.begin() + 1; pI_it != peptide_index.end(); ++pI_it)
        {
          // swap current position with random element from vector
          // swapping positions are random in range [0, current_position + 1)
          // which can be at most [0, n)
          std::iter_swap(pI_it, peptide_index.begin() + pseudoRNG((pI_it - peptide_index.begin()) + 1));
        }
      }

      // re-insert the missing K/P/R at the appropriate places
      for (IndexType::iterator it = idx.begin(); it != idx.end(); ++it)
      {
        peptide_index.insert(peptide_index.begin() + it->first, it->first);
      }

      // use the shuffled index to create the get the new peptide sequence and
      // then to place the modifications at their appropriate places (at the
      // same, shuffled AA where they were before).
      for (Size i = 0; i < peptide_index.size(); i++)
      {
        shuffled.sequence[i] = peptide.sequence[peptide_index[i]];
      }
      for (Size j = 0; j < shuffled.mods.size(); j++)
      {
        for (Size k = 0; k < peptide_index.size(); k++)
        {
          // C and N terminal mods are implicitly not shuffled because they live at positions -1 and sequence.size()
          if (boost::numeric_cast<int>(peptide_index[k]) == shuffled.mods[j].location)
          {
            shuffled.mods[j].location = boost::numeric_cast<int>(k);
            break;
          }
        }
      }

#ifdef DEBUG_MRMDECOY
      for (Size j = 0; j < shuffled.mods.size(); j++)
      {
        std::cout << " position after shuffling " << shuffled.mods[j].location << " mass difference " << shuffled.mods[j].mono_mass_delta << std::endl;
      }
#endif

      ++attempts;

      // If our attempts have failed so far, we will append two random AA to
      // the sequence and see whether we can achieve sufficient shuffling with
      // these additional AA added to the sequence.
      if (attempts % 10 == 9)
      {
        if (replace_aa_instead_append)
        {
          OpenMS::AASequence shuffled_sequence = TargetedExperimentHelper::getAASequence(shuffled);
          int res_pos = (pseudoRNG() % aa_size);
          int pep_pos = -1;
          size_t pos_trials = 0;
          while (pep_pos < 0 && pos_trials < shuffled_sequence.size())
          {
            pep_pos = (pseudoRNG() % shuffled_sequence.size());
            if (shuffled_sequence[pep_pos].isModified() || (shuffled_sequence.hasNTerminalModification() && pep_pos == 0) || (shuffled_sequence.hasNTerminalModification() && pep_pos == (int)(shuffled_sequence.size() - 1)))
            {
              pep_pos = -1;
            }
            else
            {
              if (pep_pos == 0)
              {
                shuffled_sequence = AASequence::fromString(aa[res_pos]) + shuffled_sequence.getSuffix(shuffled_sequence.size() - pep_pos - 1);
              }
              else if (pep_pos == (int)(shuffled_sequence.size() - 1))
              {
                shuffled_sequence = shuffled_sequence.getPrefix(pep_pos) + AASequence::fromString(aa[res_pos]);
              }
              else
              {
                shuffled_sequence = shuffled_sequence.getPrefix(pep_pos) + AASequence::fromString(aa[res_pos]) + shuffled_sequence.getSuffix(shuffled_sequence.size() - pep_pos - 1);
              }
            }
            ++pos_trials;
          }
          shuffled.sequence = shuffled_sequence.toUnmodifiedString();
          peptide = shuffled;
        }
        else
        {
          int pos = (pseudoRNG() % aa_size);
          peptide.sequence.append(aa[pos]);
          pos = (pseudoRNG() % aa_size);
          peptide.sequence.append(aa[pos]);
          // now make the shuffled peptide the same length as the new peptide
          shuffled = peptide;
        }
      }
    }

    return shuffled;
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::pseudoreversePeptide(
    OpenMS::TargetedExperiment::Peptide peptide)
  {
    OpenMS::TargetedExperiment::Peptide peptideorig = peptide;
    std::vector<Size> peptide_index;
    for (Size i = 0; i < peptide.sequence.size(); i++)
    {
      peptide_index.push_back(i);
    }

    peptide.sequence = peptide.sequence.substr(0, peptide.sequence.size() - 1).reverse()
                       + peptide.sequence.substr(peptide.sequence.size() - 1, 1); // pseudo-reverse
    std::reverse(peptide_index.begin(), peptide_index.end() - 1);

    for (Size j = 0; j < peptide.mods.size(); j++)
    {
      for (Size k = 0; k < peptide_index.size(); k++)
      {
        if (boost::numeric_cast<int>(peptide_index[k])  == peptide.mods[j].location)
        {
          peptide.mods[j].location = boost::numeric_cast<int>(k);
          break;
        }
      }
    }

    return peptide;
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::reversePeptide(
    OpenMS::TargetedExperiment::Peptide peptide)
  {
    OpenMS::TargetedExperiment::Peptide peptideorig = peptide;
    std::vector<Size> peptide_index;
    for (Size i = 0; i < peptide.sequence.size(); i++)
    {
      peptide_index.push_back(i);
    }

    peptide.sequence = peptide.sequence.reverse();
    std::reverse(peptide_index.begin(), peptide_index.end());

    for (Size j = 0; j < peptide.mods.size(); j++)
    {
      for (Size k = 0; k < peptide_index.size(); k++)
      {
        if (boost::numeric_cast<int>(peptide_index[k]) == peptide.mods[j].location)
        {
          peptide.mods[j].location = boost::numeric_cast<int>(k);
          break;
        }
      }
    }

    return peptide;
  }

  bool MRMDecoy::has_CNterminal_mods(const OpenMS::TargetedExperiment::Peptide& peptide)
  {
    for (Size j = 0; j < peptide.mods.size(); j++)
    {
      if (peptide.mods[j].location == -1 || peptide.mods[j].location == boost::numeric_cast<int>(peptide.sequence.size()))
      {
        return true;
      }
    }
    return false;
  }

  void MRMDecoy::generateDecoys(OpenMS::TargetedExperiment& exp, OpenMS::TargetedExperiment& dec,
                                String method, String decoy_tag, double identity_threshold, int max_attempts,
                                double mz_threshold, double mz_shift, bool exclude_similar,
                                double similarity_threshold, bool remove_CNterminal_mods, double precursor_mass_shift,
                                std::vector<String> fragment_types, std::vector<size_t> fragment_charges,
                                bool enable_specific_losses, bool enable_unspecific_losses, bool remove_unannotated,
                                int round_decPow)
  {
    MRMIonSeries mrmis;
    MRMDecoy::PeptideVectorType peptides, decoy_peptides;
    MRMDecoy::ProteinVectorType proteins, decoy_proteins;
    MRMDecoy::TransitionVectorType decoy_transitions;
    for (Size i = 0; i < exp.getProteins().size(); i++)
    {
      OpenMS::TargetedExperiment::Protein protein = exp.getProteins()[i];
      protein.id = decoy_tag + protein.id;
      proteins.push_back(protein);
    }

    std::vector<String> exclusion_peptides;
    // Go through all peptides and apply the decoy method to the sequence
    // (pseudo-reverse, reverse or shuffle). Then set the peptides and proteins of the decoy
    // experiment.
    for (Size pep_idx = 0; pep_idx < exp.getPeptides().size(); ++pep_idx)
    {
      OpenMS::TargetedExperiment::Peptide peptide = exp.getPeptides()[pep_idx];
      // continue if the peptide has C/N terminal modifications and we should exclude them
      if (remove_CNterminal_mods && MRMDecoy::has_CNterminal_mods(peptide)) {continue; }
      peptide.id = decoy_tag + peptide.id;
      OpenMS::String original_sequence = peptide.sequence;
      if (!peptide.getPeptideGroupLabel().empty())
      {
        peptide.setPeptideGroupLabel(decoy_tag + peptide.getPeptideGroupLabel());
      }

      if (method == "pseudo-reverse")
      {
        peptide = MRMDecoy::pseudoreversePeptide(peptide);
      }
      else if (method == "reverse")
      {
        peptide = MRMDecoy::reversePeptide(peptide);
      }
      else if (method == "shuffle")
      {
        peptide = MRMDecoy::shufflePeptide(peptide, identity_threshold, -1, max_attempts);
      }
      for (Size prot_idx = 0; prot_idx < peptide.protein_refs.size(); ++prot_idx)
      {
        peptide.protein_refs[prot_idx] = decoy_tag + peptide.protein_refs[prot_idx];
      }

      if (MRMDecoy::AASequenceIdentity(original_sequence, peptide.sequence) > identity_threshold)
      {
        if (!exclude_similar)
        {
          std::cout << "Target sequence: " << original_sequence << " Decoy sequence: " << peptide.sequence  << " Sequence identity: " << MRMDecoy::AASequenceIdentity(original_sequence, peptide.sequence) << " Identity threshold: " << identity_threshold << std::endl;
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "AA Sequences are too similar. Either decrease identity_threshold and increase max_attempts for the shuffle method or set flag exclude_similar.");
        }
        else
        {
          exclusion_peptides.push_back(peptide.id);
        }
      }

      peptides.push_back(peptide);
    }
    dec.setPeptides(peptides); // temporary set peptides, overwrite later again!

    // hash of the peptide reference containing all transitions
    MRMDecoy::PeptideTransitionMapType peptide_trans_map;
    for (Size i = 0; i < exp.getTransitions().size(); i++)
    {
      peptide_trans_map[exp.getTransitions()[i].getPeptideRef()].push_back(&exp.getTransitions()[i]);
    }

    Size progress = 0;
    startProgress(0, exp.getTransitions().size(), "Creating decoys");
    for (MRMDecoy::PeptideTransitionMapType::iterator pep_it = peptide_trans_map.begin();
         pep_it != peptide_trans_map.end(); ++pep_it)
    {
      String peptide_ref = pep_it->first;
      String decoy_peptide_ref = decoy_tag + pep_it->first; // see above, the decoy peptide id is computed deterministically from the target id
      const TargetedExperiment::Peptide target_peptide = exp.getPeptideByRef(peptide_ref);
      // continue if the peptide has C/N terminal modifications and we should exclude them
      if (remove_CNterminal_mods && MRMDecoy::has_CNterminal_mods(target_peptide)) {continue;}

      const TargetedExperiment::Peptide decoy_peptide = dec.getPeptideByRef(decoy_peptide_ref);
      OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);
      OpenMS::AASequence decoy_peptide_sequence = TargetedExperimentHelper::getAASequence(decoy_peptide);

      int decoy_charge = 1;
      int target_charge = 1;
      if (decoy_peptide.hasCharge()) {decoy_charge = decoy_peptide.getChargeState();}
      if (target_peptide.hasCharge()) {target_charge = target_peptide.getChargeState();}

      MRMIonSeries::IonSeries decoy_ionseries = mrmis.getIonSeries(decoy_peptide_sequence, decoy_charge, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, round_decPow);
      MRMIonSeries::IonSeries target_ionseries = mrmis.getIonSeries(target_peptide_sequence, target_charge, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, round_decPow);

      for (Size i = 0; i < pep_it->second.size(); i++)
      {
        setProgress(++progress);
        const ReactionMonitoringTransition tr = *(pep_it->second[i]);

        if (!tr.isDetectingTransition() || tr.getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
        {
          continue;
        }

        ReactionMonitoringTransition decoy_tr = tr; // copy the target transition

        decoy_tr.setNativeID(decoy_tag + tr.getNativeID());
        decoy_tr.setDecoyTransitionType(ReactionMonitoringTransition::DECOY);
        decoy_tr.setPrecursorMZ(tr.getPrecursorMZ() + precursor_mass_shift); // fix for TOPPView: Duplicate precursor MZ is not displayed.

        // determine the current annotation for the target ion and then select
        // the appropriate decoy ion for this target transition
        std::pair<String, double> targetion = mrmis.annotateIon(target_ionseries, tr.getProductMZ(), mz_threshold);
        std::pair<String, double> decoyion = mrmis.getIon(decoy_ionseries, targetion.first);

        if (method == "shift")
        {
          decoy_tr.setProductMZ(decoyion.second + mz_shift);
        }
        else
        {
          decoy_tr.setProductMZ(decoyion.second);
        }
        decoy_tr.setPeptideRef(decoy_tag + tr.getPeptideRef());

        if (decoyion.second > 0)
        {
          if (similarity_threshold >= 0)
          {
            if (std::fabs(tr.getProductMZ() - decoy_tr.getProductMZ()) < similarity_threshold)
            {
              if (!exclude_similar)
              {
                throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Fragment ions are too similar. Either decrease similarity_threshold or set flag exclude_similar.");
              }
              else
              {
                exclusion_peptides.push_back(decoy_tr.getPeptideRef());
              } 
           }
          }
          decoy_transitions.push_back(decoy_tr);
        }
        else
        {
          if (remove_unannotated)
          {
            exclusion_peptides.push_back(decoy_tr.getPeptideRef());
          }
          else
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decoy fragment ion for target fragment ion " + String(targetion.first) + " of peptide " + target_peptide_sequence.toString() + " with precursor charge " + String(target_peptide.getChargeState()) + " could not be mapped. Please check whether it is a valid ion and enable losses or removal of terminal modifications if necessary. Skipping of unannotated target assays is available as last resort.");
          }
        }
      } // end loop over transitions
    } // end loop over peptides
    endProgress();

    MRMDecoy::TransitionVectorType filtered_decoy_transitions;
    for (MRMDecoy::TransitionVectorType::iterator tr_it = decoy_transitions.begin(); tr_it != decoy_transitions.end(); ++tr_it)
    {
      if (std::find(exclusion_peptides.begin(), exclusion_peptides.end(), tr_it->getPeptideRef()) == exclusion_peptides.end())
      {
        filtered_decoy_transitions.push_back(*tr_it);
      }
    }
    dec.setTransitions(filtered_decoy_transitions);

    std::vector<String> protein_ids;
    for (Size i = 0; i < peptides.size(); ++i)
    {
      TargetedExperiment::Peptide peptide = peptides[i];

      // Check if peptide has any transitions left
      if (std::find(exclusion_peptides.begin(), exclusion_peptides.end(), peptide.id) == exclusion_peptides.end())
      {
        decoy_peptides.push_back(peptide);
        for (Size j = 0; j < peptide.protein_refs.size(); ++j)
        {
          protein_ids.push_back(peptide.protein_refs[j]);
        }
      }
      else
      {
        LOG_DEBUG << "[peptide] Skipping " << peptide.id << std::endl;
      }
    }

    for (Size i = 0; i < proteins.size(); ++i)
    {
      OpenMS::TargetedExperiment::Protein protein = proteins[i];

      // Check if protein has any peptides left
      if (find(protein_ids.begin(), protein_ids.end(), protein.id) != protein_ids.end())
      {
        decoy_proteins.push_back(protein);
      }
      else
      {
        LOG_DEBUG << "[protein] Skipping " << protein.id << std::endl;
      }
    }

    dec.setPeptides(decoy_peptides);
    dec.setProteins(decoy_proteins);

  }

}

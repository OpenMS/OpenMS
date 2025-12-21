// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/unordered_map.hpp>
#include <unordered_set>

namespace OpenMS
{

  MRMDecoy::MRMDecoy() :
    DefaultParamHandler("MRMDecoy"),
    ProgressLogger()
  {
    defaults_.setValue("non_shuffle_pattern", "KRP", "Residues to not shuffle (keep at a constant position when shuffling). Default is 'KPR' to not shuffle lysine, arginine and proline.");

    defaults_.setValue("keepPeptideNTerm", "true", "Whether to keep peptide N terminus constant when shuffling / reversing.", {"advanced"});
    defaults_.setValidStrings("keepPeptideNTerm", {"true","false"});

    defaults_.setValue("keepPeptideCTerm", "true", "Whether to keep peptide C terminus constant when shuffling / reversing.", {"advanced"});
    defaults_.setValidStrings("keepPeptideCTerm", {"true","false"});

    // write defaults into Param object param_
    defaultsToParam_();
  }

  void MRMDecoy::updateMembers_()
  {
    keep_const_pattern_ = param_.getValue("non_shuffle_pattern").toString();
    keepN_ = param_.getValue("keepPeptideNTerm").toBool();
    keepC_ = param_.getValue("keepPeptideCTerm").toBool();
  }

  MRMDecoy::IndexType MRMDecoy::findFixedResidues(const std::string& sequence,
      bool keepN, bool keepC, const OpenMS::String& keep_const_pattern)
  {
    // also blocks both N- and C-terminus from shuffling if required
    MRMDecoy::IndexType idx;
    for (size_t i = 0; i < sequence.size(); i++)
    {
      if ( (keepN && i == 0) || (keepC && i + 1 == sequence.size()))
      {
        idx.push_back(i);
        continue;
      }

      for (size_t j = 0; j < keep_const_pattern.size(); j++)
      {
        if (sequence[i] == keep_const_pattern[j])
        {
          idx.push_back(i);
        }
      }
    }
    return idx;
  }

  MRMDecoy::IndexType MRMDecoy::findFixedResidues_(const std::string& sequence) const
  {
    return MRMDecoy::findFixedResidues(sequence, false, false, keep_const_pattern_);
  }

  MRMDecoy::IndexType MRMDecoy::findFixedAndTermResidues_(const std::string& sequence) const
  {
    return MRMDecoy::findFixedResidues(sequence, keepN_, keepC_, keep_const_pattern_);
  }

  float MRMDecoy::AASequenceIdentity(const String& sequence, const String& decoy) const
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
    OpenMS::TargetedExperiment::Peptide peptide, const double identity_threshold, int seed,
    const int max_attempts) const
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

    //TODO use Math::RandomShuffler (=same approach) and put this rather general method somewhere more prominent.
    boost::mt19937 generator(seed);
    boost::uniform_int<> uni_dist;
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > pseudoRNG(generator, uni_dist);

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
      // Block tryptic residues and N-/C-terminus from shuffling
      MRMDecoy::IndexType idx = findFixedAndTermResidues_(peptide.sequence);

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
        peptide_index.erase(peptide_index.begin() + *it);
      }

      // shuffle the peptide index (without the K/P/R which we leave in place)
      // one could also use std::random_shuffle here but then the code becomes
      // not testable since the implementation of std::random_shuffle differs
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
        peptide_index.insert(peptide_index.begin() + *it, *it);
      }

      // use the shuffled index to create the new peptide sequence and
      // then to place the modifications at their appropriate places (make sure
      // that the modifications are placed with their initial amino acids).
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

      // If our attempts have failed so far, we will mutate a random AA of
      // the sequence and see whether we can achieve sufficient shuffling with
      // the new sequence.
      if (attempts % 10 == 9)
      {
        OpenMS::AASequence shuffled_sequence = TargetedExperimentHelper::getAASequence(shuffled);
        int res_pos = (pseudoRNG() % aa_size);
        int pep_pos = -1;
        size_t pos_trials = 0;
        while (pep_pos < 0 && pos_trials < shuffled_sequence.size())
        {
          // select position to mutate (and ensure we are not changing N/C terminus or any modified position doing it)
          pep_pos = (pseudoRNG() % shuffled_sequence.size());
          if (shuffled_sequence[pep_pos].isModified() || (pep_pos == 0) || (pep_pos == (int)(shuffled_sequence.size() - 1)))
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
    }

    return shuffled;
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::reversePeptide(
      const OpenMS::TargetedExperiment::Peptide& peptide, const bool keepN, const bool keepC, 
      const String& const_pattern)
  {
    OpenMS::TargetedExperiment::Peptide reversed = peptide;
    // Block tryptic residues and N-/C-terminus from shuffling
    MRMDecoy::IndexType idx = MRMDecoy::findFixedResidues(peptide.sequence, keepN, keepC, const_pattern);

    std::vector<Size> peptide_index;
    for (Size i = 0; i < peptide.sequence.size(); i++)
    {
      peptide_index.push_back(i);
    }

    // we erase the indices where K/P/R are (from the back / in reverse order
    // to not delete indices we access later)
    for (IndexType::reverse_iterator it = idx.rbegin(); it != idx.rend(); ++it)
    {
      peptide_index.erase(peptide_index.begin() + *it);
    }

    // reverse the peptide index
    std::reverse(peptide_index.begin(), peptide_index.end());

    // re-insert the missing K/P/R at the appropriate places
    for (IndexType::iterator it = idx.begin(); it != idx.end(); ++it)
    {
      peptide_index.insert(peptide_index.begin() + *it, *it);
    }

    // use the reversed index to create the new peptide sequence and
    // then to place the modifications at their appropriate places (make sure
    // that the modifications are placed with their initial amino acids).
    for (Size i = 0; i < peptide_index.size(); i++)
    {
      reversed.sequence[i] = peptide.sequence[peptide_index[i]];
    }
    for (Size j = 0; j < reversed.mods.size(); j++)
    {
      for (Size k = 0; k < peptide_index.size(); k++)
      {
        // C and N terminal mods are implicitly not reversed because they live at positions -1 and sequence.size()
        if (boost::numeric_cast<int>(peptide_index[k]) == reversed.mods[j].location)
        {
          reversed.mods[j].location = boost::numeric_cast<int>(k);
          break;
        }
      }
    }
    return reversed;
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::pseudoreversePeptide_(
    const OpenMS::TargetedExperiment::Peptide& peptide) const
  {
    return MRMDecoy::reversePeptide(peptide, false, true);
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::reversePeptide_(
    const OpenMS::TargetedExperiment::Peptide& peptide) const
  {
    return MRMDecoy::reversePeptide(peptide, false, false);
  }

  void MRMDecoy::switchKR(OpenMS::TargetedExperiment::Peptide& peptide) const
  {
    static std::string aa[] =
    {
      "A", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "M", "F", "S", "T", "W",
      "Y", "V"
    };
    int aa_size = 17;

    static boost::mt19937 generator(42);
    static boost::uniform_int<> uni_dist;
    static boost::variate_generator<boost::mt19937&, boost::uniform_int<> > pseudoRNG(generator, uni_dist);

    Size lastAA = peptide.sequence.size() - 1;
    if (peptide.sequence[lastAA] == 'K')
    {
      peptide.sequence[lastAA] = 'R';
    }
    else if (peptide.sequence[lastAA] == 'R')
    {
      peptide.sequence[lastAA] = 'K';
    }
    else
    {
      // randomize
      int res_pos = (pseudoRNG() % aa_size);
      peptide.sequence[lastAA] = (char)aa[res_pos][0];
    }
  }

  bool MRMDecoy::hasCNterminalMods_(const OpenMS::TargetedExperiment::Peptide& peptide, bool checkCterminalAA) const
  {
    for (Size j = 0; j < peptide.mods.size(); j++)
    {
      if (peptide.mods[j].location == -1 || peptide.mods[j].location == (int)peptide.sequence.size())
      {
        return true;
      }
      if (checkCterminalAA && peptide.mods[j].location == (int)peptide.sequence.size() - 1)
      {
        return true;
      }
    }
    return false;
  }

  String MRMDecoy::getModifiedPeptideSequence_(const OpenMS::TargetedExperiment::Peptide& pep) const
  {
    String full_peptide_name;
    for (int loc = -1; loc <= (int)pep.sequence.size(); loc++)
    {
      if (loc > -1 && loc < (int)pep.sequence.size())
      {
        full_peptide_name += pep.sequence[loc];
      }
      // C-terminal and N-terminal modifications may be at positions -1 or pep.sequence
      for (Size modloc = 0; modloc < pep.mods.size(); modloc++)
      {
        if (pep.mods[modloc].location == loc)
        {
          full_peptide_name += "(UniMod:" + String(pep.mods[modloc].unimod_id) + ")";
        }
      }
    }
    return full_peptide_name;
  }
  
  void MRMDecoy::generateDecoys(const OpenMS::TargetedExperiment& exp, OpenMS::TargetedExperiment& dec,
                                const String& method, const double aim_decoy_fraction, const bool do_switchKR,
                                const String& decoy_tag, const int max_attempts, const double identity_threshold,
                                const double precursor_mz_shift, const double product_mz_shift, const double product_mz_threshold,
                                const std::vector<String>& fragment_types, const std::vector<size_t>& fragment_charges,
                                const bool enable_specific_losses, const bool enable_unspecific_losses, const int round_decPow) const
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

    srand(time(nullptr));
    std::vector<size_t> item_list, selection_list;
    item_list.reserve(exp.getPeptides().size());
    for (Size k = 0; k < exp.getPeptides().size(); k++) {item_list.push_back(k);}

    if (aim_decoy_fraction > 1.0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decoy fraction needs to be less than one (values larger than one currently not supported).");
    }
    else if (aim_decoy_fraction < 1.0)
    {
      Math::RandomShuffler shuffler;
      shuffler.portable_random_shuffle(item_list.begin(), item_list.end());
      selection_list.reserve(aim_decoy_fraction * exp.getPeptides().size());
      Size k = 0;
      while (selection_list.size() < aim_decoy_fraction * exp.getPeptides().size())
      {
        selection_list.push_back( item_list[ k++ % item_list.size() ]);
      }
    }
    else
    {
      selection_list = item_list;
    }

    // Create an unordered_map of all precursors (modified peptide sequence + charge) to their IDs, this will be used to search if a decoy sequence is also a target sequence
    std::unordered_map<std::string, std::string> allPeptideSequences;
    for (const auto& pep_idx: selection_list)
    {
      OpenMS::TargetedExperiment::Peptide peptide = exp.getPeptides()[pep_idx];

      // create a modified peptide sequence string
      allPeptideSequences[MRMDecoy::getModifiedPeptideSequence_(peptide) + String(peptide.getChargeState()) ] = peptide.id;
    }

    std::unordered_set<String> exclusion_peptides;
    // Go through all peptides and apply the decoy method to the sequence
    // (pseudo-reverse, reverse or shuffle). Then set the peptides and proteins of the decoy
    // experiment.
    Size progress = 0;
    startProgress(0, selection_list.size(), "Generating decoy peptides");
    for (const auto& pep_idx : selection_list)
    {
      setProgress(++progress);

      OpenMS::TargetedExperiment::Peptide peptide = exp.getPeptides()[pep_idx];

      peptide.id = decoy_tag + peptide.id;

      if (!peptide.getPeptideGroupLabel().empty())
      {
        peptide.setPeptideGroupLabel(decoy_tag + peptide.getPeptideGroupLabel());
      }

      if (method == "pseudo-reverse")
      {
        // exclude peptide if it has C/N terminal modifications because we can't do a (partial) reverse
        if (MRMDecoy::hasCNterminalMods_(peptide, do_switchKR))
        {
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " due to C/N-terminal modifications" << std::endl;
          exclusion_peptides.insert(peptide.id);
        }
        else
        {
          peptide = MRMDecoy::pseudoreversePeptide_(peptide);
          if (do_switchKR) { switchKR(peptide); }
        }
      }
      else if (method == "reverse")
      {
        // exclude peptide if it has C/N terminal modifications because we can't do a (partial) reverse
        if (MRMDecoy::hasCNterminalMods_(peptide, false))
        {
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " due to C/N-terminal modifications" << std::endl;
          exclusion_peptides.insert(peptide.id);
        }
        else
        {
          peptide = MRMDecoy::reversePeptide_(peptide);
        }
      }
      else if (method == "shuffle")
      {
        peptide = MRMDecoy::shufflePeptide(peptide, identity_threshold, -1, max_attempts);
        if (do_switchKR && MRMDecoy::hasCNterminalMods_(peptide, do_switchKR))
        {
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " due to C/N-terminal modifications" << std::endl;
          exclusion_peptides.insert(peptide.id);
        }
        else if (do_switchKR) { switchKR(peptide); }
      }

      // Check that the decoy precursor does not happen to be a target precursor AND is not already present
      if (allPeptideSequences.find(peptide.sequence + String((peptide.getChargeState()))) != allPeptideSequences.end())
      {
        OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " since decoy peptide is also a target peptide or this decoy peptide is already present" << std::endl;
        exclusion_peptides.insert(peptide.id);
      }
      else
      {
        // Since this decoy will be added, add it to the precursor map so that the same decoy is not added twice
        OPENMS_LOG_DEBUG << "[peptide] adding " << peptide.id << " to master list of peptides " << std::endl;
        allPeptideSequences[MRMDecoy::getModifiedPeptideSequence_(peptide) + String(peptide.getChargeState())] = peptide.id;
      }

      for (Size prot_idx = 0; prot_idx < peptide.protein_refs.size(); ++prot_idx)
      {
        peptide.protein_refs[prot_idx] = decoy_tag + peptide.protein_refs[prot_idx];
      }

      peptides.push_back(peptide);
    }
    endProgress();
    dec.setPeptides(peptides); // temporary set peptides, overwrite later again!

    // hash of the peptide reference containing all transitions
    MRMDecoy::PeptideTransitionMapType peptide_trans_map;
    for (Size i = 0; i < exp.getTransitions().size(); i++)
    {
      peptide_trans_map[exp.getTransitions()[i].getPeptideRef()].push_back(&exp.getTransitions()[i]);
    }

    progress = 0;
    startProgress(0, peptide_trans_map.size(), "Generating decoy transitions");
    for (MRMDecoy::PeptideTransitionMapType::iterator pep_it = peptide_trans_map.begin();
         pep_it != peptide_trans_map.end(); ++pep_it)
    {
      setProgress(++progress);

      String peptide_ref = pep_it->first;
      String decoy_peptide_ref = decoy_tag + pep_it->first; // see above, the decoy peptide id is computed deterministically from the target id
      if (!dec.hasPeptide(decoy_peptide_ref)) { continue; }
      const TargetedExperiment::Peptide& target_peptide = exp.getPeptideByRef(peptide_ref);

      const TargetedExperiment::Peptide& decoy_peptide = dec.getPeptideByRef(decoy_peptide_ref);
      OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);
      OpenMS::AASequence decoy_peptide_sequence = TargetedExperimentHelper::getAASequence(decoy_peptide);

      int decoy_charge = 1;
      int target_charge = 1;
      if (decoy_peptide.hasCharge()) { decoy_charge = decoy_peptide.getChargeState(); }
      if (target_peptide.hasCharge()) { target_charge = target_peptide.getChargeState(); }

      MRMIonSeries::IonSeries decoy_ionseries = mrmis.getIonSeries(decoy_peptide_sequence, decoy_charge,
                                                                   fragment_types, fragment_charges, enable_specific_losses,
                                                                   enable_unspecific_losses, round_decPow);
      MRMIonSeries::IonSeries target_ionseries = mrmis.getIonSeries(target_peptide_sequence, target_charge,
                                                                    fragment_types, fragment_charges, enable_specific_losses,
                                                                    enable_unspecific_losses, round_decPow);

      // Compute (new) decoy precursor m/z based on the K/R replacement and the AA changes in the shuffle algorithm
      double decoy_precursor_mz = decoy_peptide_sequence.getMZ(decoy_charge);
      decoy_precursor_mz += precursor_mz_shift; // fix for TOPPView: Duplicate precursor MZ is not displayed.

      for (Size i = 0; i < pep_it->second.size(); i++)
      {
        const ReactionMonitoringTransition tr = *(pep_it->second[i]);

        if (!tr.isDetectingTransition() || tr.getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
        {
          continue;
        }

        ReactionMonitoringTransition decoy_tr = tr; // copy the target transition

        decoy_tr.setNativeID(decoy_tag + tr.getNativeID());
        decoy_tr.setDecoyTransitionType(ReactionMonitoringTransition::DECOY);
        decoy_tr.setPrecursorMZ(decoy_precursor_mz);

        // determine the current annotation for the target ion and then select
        // the appropriate decoy ion for this target transition
        std::pair<String, double> targetion = mrmis.annotateIon(target_ionseries, tr.getProductMZ(), product_mz_threshold);
        std::pair<String, double> decoyion = mrmis.getIon(decoy_ionseries, targetion.first);

        if (method == "shift")
        {
          decoy_tr.setProductMZ(decoyion.second + product_mz_shift);
        }
        else
        {
          decoy_tr.setProductMZ(decoyion.second);
        }
        decoy_tr.setPeptideRef(decoy_tag + tr.getPeptideRef());

        if (decoyion.second > 0)
        {
          decoy_transitions.push_back(decoy_tr);
        }
        else
        {
          // transition could not be annotated, remove whole peptide
          exclusion_peptides.insert(decoy_tr.getPeptideRef());
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << decoy_tr.getPeptideRef() << " due to missing annotation" << std::endl;
        }
      } // end loop over transitions

    } // end loop over peptides
    endProgress();

    decoy_transitions.erase(std::remove_if(
                            decoy_transitions.begin(), decoy_transitions.end(),
                            [&exclusion_peptides](const OpenMS::ReactionMonitoringTransition& tr)
    {
      return exclusion_peptides.find(tr.getPeptideRef()) != exclusion_peptides.end();
    }), decoy_transitions.end());
    dec.setTransitions(std::move(decoy_transitions));

    std::unordered_set<String> protein_ids;
    decoy_peptides.reserve(peptides.size());
    for (const auto& peptide : peptides)
    {
      // Check if peptide has any transitions left
      if (exclusion_peptides.find(peptide.id) == exclusion_peptides.end())
      {
        decoy_peptides.push_back(peptide);
        for (Size j = 0; j < peptide.protein_refs.size(); ++j)
        {
          protein_ids.insert(peptide.protein_refs[j]);
        }
      }
      else
      {
        OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " due to missing transitions" << std::endl;
      }
    }

    decoy_proteins.reserve(proteins.size());
    for (const auto& protein : proteins)
    {
      // Check if protein has any peptides left
      if (protein_ids.find(protein.id) != protein_ids.end())
      {
        decoy_proteins.push_back(protein);
      }
      else
      {
        OPENMS_LOG_DEBUG << "[protein] Skipping " << protein.id << " due to missing peptides" << std::endl;
      }
    }

    dec.setPeptides(std::move(decoy_peptides));
    dec.setProteins(std::move(decoy_proteins));

  }

  void MRMDecoy::generateDecoysLight(const OpenSwath::LightTargetedExperiment& exp,
                                     OpenSwath::LightTargetedExperiment& dec,
                                     const String& method,
                                     const double aim_decoy_fraction,
                                     const bool do_switchKR,
                                     const String& decoy_tag,
                                     const int max_attempts,
                                     const double identity_threshold,
                                     const double precursor_mz_shift,
                                     const double product_mz_shift,
                                     const double product_mz_threshold,
                                     const std::vector<String>& fragment_types,
                                     const std::vector<size_t>& fragment_charges,
                                     const bool enable_specific_losses,
                                     const bool enable_unspecific_losses,
                                     const int round_decPow) const
  {
    MRMIonSeries mrmis;
    std::vector<OpenSwath::LightCompound> decoy_compounds;
    std::vector<OpenSwath::LightProtein> decoy_proteins;
    std::vector<OpenSwath::LightTransition> decoy_transitions;

    // Create decoy proteins
    for (const auto& protein : exp.proteins)
    {
      OpenSwath::LightProtein decoy_protein = protein;
      decoy_protein.id = decoy_tag + protein.id;
      decoy_proteins.push_back(decoy_protein);
    }

    srand(time(nullptr));
    std::vector<size_t> item_list, selection_list;
    item_list.reserve(exp.compounds.size());
    for (Size k = 0; k < exp.compounds.size(); k++) { item_list.push_back(k); }

    if (aim_decoy_fraction > 1.0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Decoy fraction needs to be less than one");
    }
    else if (aim_decoy_fraction < 1.0)
    {
      Math::RandomShuffler shuffler;
      shuffler.portable_random_shuffle(item_list.begin(), item_list.end());
      selection_list.reserve(aim_decoy_fraction * exp.compounds.size());
      Size k = 0;
      while (selection_list.size() < aim_decoy_fraction * exp.compounds.size())
      {
        selection_list.push_back(item_list[k++ % item_list.size()]);
      }
    }
    else
    {
      selection_list = item_list;
    }

    // Create map of all peptide sequences to detect duplicates
    std::unordered_map<std::string, std::string> allPeptideSequences;
    for (const auto& idx : selection_list)
    {
      const auto& compound = exp.compounds[idx];
      if (compound.isPeptide())
      {
        allPeptideSequences[compound.sequence + std::to_string(compound.charge)] = compound.id;
      }
    }

    std::unordered_set<std::string> exclusion_peptides;

    // Build compound map
    std::map<std::string, const OpenSwath::LightCompound*> compound_map;
    for (const auto& compound : exp.compounds)
    {
      compound_map[compound.id] = &compound;
    }

    // Go through all compounds and generate decoys
    Size progress = 0;
    startProgress(0, selection_list.size(), "Generating decoy peptides (Light)");

    for (const auto& idx : selection_list)
    {
      setProgress(++progress);

      OpenSwath::LightCompound decoy_compound = exp.compounds[idx];
      const std::string original_id = decoy_compound.id;
      decoy_compound.id = decoy_tag + decoy_compound.id;

      // Update protein refs
      for (auto& prot_ref : decoy_compound.protein_refs)
      {
        prot_ref = decoy_tag + prot_ref;
      }

      // Update peptide_group_label
      if (!decoy_compound.peptide_group_label.empty())
      {
        decoy_compound.peptide_group_label = decoy_tag + decoy_compound.peptide_group_label;
      }

      if (!decoy_compound.isPeptide())
      {
        // For metabolites, just copy with decoy tag
        decoy_compounds.push_back(decoy_compound);
        continue;
      }

      // Parse sequence to AASequence
      OpenMS::AASequence original_sequence;
      try
      {
        original_sequence = AASequence::fromString(decoy_compound.sequence);
      }
      catch (Exception::InvalidValue&)
      {
        OPENMS_LOG_DEBUG << "[peptide] Skipping " << decoy_compound.id << " - cannot parse sequence" << std::endl;
        exclusion_peptides.insert(decoy_compound.id);
        continue;
      }

      // Create a temporary heavy Peptide for decoy generation
      OpenMS::TargetedExperiment::Peptide temp_peptide;
      temp_peptide.sequence = original_sequence.toUnmodifiedString();
      temp_peptide.id = original_id;

      // Copy modifications
      for (const auto& mod : decoy_compound.modifications)
      {
        OpenMS::TargetedExperiment::Peptide::Modification temp_mod;
        temp_mod.location = mod.location;
        temp_mod.unimod_id = mod.unimod_id;
        temp_peptide.mods.push_back(temp_mod);
      }

      // Apply decoy method
      OpenMS::TargetedExperiment::Peptide decoy_temp_peptide;
      if (method == "pseudo-reverse")
      {
        if (hasCNterminalMods_(temp_peptide, do_switchKR))
        {
          exclusion_peptides.insert(decoy_compound.id);
          continue;
        }
        decoy_temp_peptide = pseudoreversePeptide_(temp_peptide);
        if (do_switchKR) { switchKR(decoy_temp_peptide); }
      }
      else if (method == "reverse")
      {
        if (hasCNterminalMods_(temp_peptide, false))
        {
          exclusion_peptides.insert(decoy_compound.id);
          continue;
        }
        decoy_temp_peptide = reversePeptide_(temp_peptide);
      }
      else if (method == "shuffle")
      {
        decoy_temp_peptide = shufflePeptide(temp_peptide, identity_threshold, -1, max_attempts);
        if (do_switchKR && hasCNterminalMods_(temp_peptide, do_switchKR))
        {
          exclusion_peptides.insert(decoy_compound.id);
          continue;
        }
        else if (do_switchKR)
        {
          switchKR(decoy_temp_peptide);
        }
      }
      else
      {
        decoy_temp_peptide = temp_peptide;
      }

      // Check for duplicates
      std::string decoy_key = std::string(decoy_temp_peptide.sequence) + std::to_string(decoy_compound.charge);
      if (allPeptideSequences.find(decoy_key) != allPeptideSequences.end())
      {
        exclusion_peptides.insert(decoy_compound.id);
        continue;
      }
      allPeptideSequences[decoy_key] = decoy_compound.id;

      // Update Light compound sequence from decoy
      decoy_compound.sequence = getModifiedPeptideSequence_(decoy_temp_peptide);

      // Update modifications
      decoy_compound.modifications.clear();
      for (const auto& mod : decoy_temp_peptide.mods)
      {
        OpenSwath::LightModification light_mod;
        light_mod.location = mod.location;
        light_mod.unimod_id = mod.unimod_id;
        decoy_compound.modifications.push_back(light_mod);
      }

      decoy_compounds.push_back(decoy_compound);
    }
    endProgress();

    // Build decoy compound map
    std::map<std::string, const OpenSwath::LightCompound*> decoy_compound_map;
    for (const auto& compound : decoy_compounds)
    {
      decoy_compound_map[compound.id] = &compound;
    }

    // Group transitions by peptide_ref
    std::map<std::string, std::vector<const OpenSwath::LightTransition*>> peptide_trans_map;
    for (const auto& tr : exp.transitions)
    {
      peptide_trans_map[tr.peptide_ref].push_back(&tr);
    }

    progress = 0;
    startProgress(0, peptide_trans_map.size(), "Generating decoy transitions (Light)");

    for (const auto& pep_it : peptide_trans_map)
    {
      setProgress(++progress);

      const std::string& peptide_ref = pep_it.first;
      std::string decoy_peptide_ref = decoy_tag + peptide_ref;

      // Check if we have this decoy compound
      auto decoy_comp_it = decoy_compound_map.find(decoy_peptide_ref);
      if (decoy_comp_it == decoy_compound_map.end())
      {
        continue;
      }

      auto target_comp_it = compound_map.find(peptide_ref);
      if (target_comp_it == compound_map.end())
      {
        continue;
      }

      const OpenSwath::LightCompound* target_compound = target_comp_it->second;
      const OpenSwath::LightCompound* decoy_compound = decoy_comp_it->second;

      if (!target_compound->isPeptide())
      {
        // For metabolites, just copy transitions with decoy tag
        for (const auto* tr : pep_it.second)
        {
          OpenSwath::LightTransition decoy_tr = *tr;
          decoy_tr.transition_name = decoy_tag + tr->transition_name;
          decoy_tr.peptide_ref = decoy_peptide_ref;
          decoy_tr.decoy = true;
          decoy_transitions.push_back(decoy_tr);
        }
        continue;
      }

      OpenMS::AASequence target_sequence, decoy_sequence;
      try
      {
        target_sequence = AASequence::fromString(target_compound->sequence);
        decoy_sequence = AASequence::fromString(decoy_compound->sequence);
      }
      catch (Exception::InvalidValue&)
      {
        continue;
      }

      int target_charge = target_compound->charge > 0 ? target_compound->charge : 1;
      int decoy_charge = decoy_compound->charge > 0 ? decoy_compound->charge : 1;

      MRMIonSeries::IonSeries target_ionseries = mrmis.getIonSeries(
          target_sequence, target_charge, fragment_types, fragment_charges,
          enable_specific_losses, enable_unspecific_losses, round_decPow);

      MRMIonSeries::IonSeries decoy_ionseries = mrmis.getIonSeries(
          decoy_sequence, decoy_charge, fragment_types, fragment_charges,
          enable_specific_losses, enable_unspecific_losses, round_decPow);

      double decoy_precursor_mz = decoy_sequence.getMZ(decoy_charge) + precursor_mz_shift;

      for (const auto* tr : pep_it.second)
      {
        if (!tr->detecting_transition || tr->decoy)
        {
          continue;
        }

        OpenSwath::LightTransition decoy_tr = *tr;
        decoy_tr.transition_name = decoy_tag + tr->transition_name;
        decoy_tr.peptide_ref = decoy_peptide_ref;
        decoy_tr.decoy = true;
        decoy_tr.precursor_mz = decoy_precursor_mz;

        // Annotate and get decoy fragment m/z
        std::pair<String, double> targetion = mrmis.annotateIon(target_ionseries, tr->product_mz, product_mz_threshold);
        std::pair<String, double> decoyion = mrmis.getIon(decoy_ionseries, targetion.first);

        if (method == "shift")
        {
          decoy_tr.product_mz = decoyion.second + product_mz_shift;
        }
        else
        {
          decoy_tr.product_mz = decoyion.second;
        }

        if (decoyion.second > 0)
        {
          // Update fragment annotation
          if (!targetion.first.empty())
          {
            decoy_tr.fragment_type = targetion.first.substr(0, 1);
            std::string num_str;
            for (size_t i = 1; i < targetion.first.size() && std::isdigit(targetion.first[i]); ++i)
            {
              num_str += targetion.first[i];
            }
            if (!num_str.empty())
            {
              decoy_tr.fragment_nr = std::stoi(num_str);
            }
            decoy_tr.annotation = targetion.first;
          }
          decoy_transitions.push_back(decoy_tr);
        }
        else
        {
          exclusion_peptides.insert(decoy_peptide_ref);
        }
      }
    }
    endProgress();

    // Filter out excluded peptides from transitions
    decoy_transitions.erase(
        std::remove_if(decoy_transitions.begin(), decoy_transitions.end(),
                       [&exclusion_peptides](const OpenSwath::LightTransition& tr) {
                         return exclusion_peptides.find(tr.peptide_ref) != exclusion_peptides.end();
                       }),
        decoy_transitions.end());

    // Filter compounds
    std::vector<OpenSwath::LightCompound> filtered_compounds;
    std::unordered_set<std::string> protein_ids;
    for (const auto& compound : decoy_compounds)
    {
      if (exclusion_peptides.find(compound.id) == exclusion_peptides.end())
      {
        filtered_compounds.push_back(compound);
        for (const auto& prot_ref : compound.protein_refs)
        {
          protein_ids.insert(prot_ref);
        }
      }
    }

    // Filter proteins
    std::vector<OpenSwath::LightProtein> filtered_proteins;
    for (const auto& protein : decoy_proteins)
    {
      if (protein_ids.find(protein.id) != protein_ids.end())
      {
        filtered_proteins.push_back(protein);
      }
    }

    dec.transitions = std::move(decoy_transitions);
    dec.compounds = std::move(filtered_compounds);
    dec.proteins = std::move(filtered_proteins);
  }

}

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

  //
  // Helper functions for converting between heavy and light modification formats
  //
  namespace detail
  {
    /// Convert heavy Peptide modifications to light format (location + unimod_id only)
    inline std::vector<OpenSwath::LightModification> toLightMods(
        const std::vector<TargetedExperiment::Peptide::Modification>& heavy_mods)
    {
      std::vector<OpenSwath::LightModification> light_mods;
      light_mods.reserve(heavy_mods.size());
      for (const auto& mod : heavy_mods)
      {
        OpenSwath::LightModification lm;
        lm.location = mod.location;
        lm.unimod_id = mod.unimod_id;
        light_mods.push_back(lm);
      }
      return light_mods;
    }

    /// Apply location updates from light modifications back to heavy modifications
    inline void applyLocationUpdates(
        std::vector<TargetedExperiment::Peptide::Modification>& heavy_mods,
        const std::vector<OpenSwath::LightModification>& light_mods)
    {
      OPENMS_PRECONDITION(heavy_mods.size() == light_mods.size(), "Modification count mismatch");
      for (Size i = 0; i < heavy_mods.size(); ++i)
      {
        heavy_mods[i].location = light_mods[i].location;
      }
    }
  } // namespace detail

  OpenMS::TargetedExperiment::Peptide MRMDecoy::shufflePeptide(
    OpenMS::TargetedExperiment::Peptide peptide, const double identity_threshold, int seed,
    const int max_attempts) const
  {
#ifdef DEBUG_MRMDECOY
    std::cout << " shuffle peptide " << peptide.sequence << '\n';
    seed = 41;
#endif
    // Delegate to light implementation
    auto light_mods = detail::toLightMods(peptide.mods);
    auto [new_sequence, new_mods] = shufflePeptideLight(
        peptide.sequence, light_mods, identity_threshold, seed, max_attempts);

    peptide.sequence = new_sequence;
    detail::applyLocationUpdates(peptide.mods, new_mods);

#ifdef DEBUG_MRMDECOY
    for (Size j = 0; j < peptide.mods.size(); j++)
    {
      std::cout << " position after shuffling " << peptide.mods[j].location << " mass difference " << peptide.mods[j].mono_mass_delta << '\n';
    }
#endif

    return peptide;
  }

  OpenMS::TargetedExperiment::Peptide MRMDecoy::reversePeptide(
      const OpenMS::TargetedExperiment::Peptide& peptide, const bool keepN, const bool keepC,
      const String& const_pattern)
  {
    // Delegate to light implementation
    auto light_mods = detail::toLightMods(peptide.mods);
    auto [new_sequence, new_mods] = reversePeptideLight(
        peptide.sequence, light_mods, keepN, keepC, const_pattern);

    OpenMS::TargetedExperiment::Peptide reversed = peptide;
    reversed.sequence = new_sequence;
    detail::applyLocationUpdates(reversed.mods, new_mods);
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
    // Delegate to light implementation (operates directly on sequence string)
    switchKRLight(peptide.sequence);
  }

  bool MRMDecoy::hasCNterminalMods_(const OpenMS::TargetedExperiment::Peptide& peptide, bool checkCterminalAA) const
  {
    // Delegate to light implementation
    auto light_mods = detail::toLightMods(peptide.mods);
    return hasCNterminalModsLight_(light_mods, peptide.sequence.size(), checkCterminalAA);
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

  //
  // Light versions of decoy generation methods
  // These operate directly on strings and LightModification vectors for memory efficiency
  //

  std::pair<std::string, std::vector<OpenSwath::LightModification>> MRMDecoy::reversePeptideLight(
      const std::string& sequence,
      const std::vector<OpenSwath::LightModification>& modifications,
      const bool keepN,
      const bool keepC,
      const String& const_pattern)
  {
    std::string reversed = sequence;
    std::vector<OpenSwath::LightModification> reversed_mods = modifications;

    // Block tryptic residues and N-/C-terminus from reversing
    MRMDecoy::IndexType idx = MRMDecoy::findFixedResidues(sequence, keepN, keepC, const_pattern);

    std::vector<Size> peptide_index;
    for (Size i = 0; i < sequence.size(); i++)
    {
      peptide_index.push_back(i);
    }

    // Erase the indices where K/P/R are (from the back to preserve indices)
    for (IndexType::reverse_iterator it = idx.rbegin(); it != idx.rend(); ++it)
    {
      peptide_index.erase(peptide_index.begin() + *it);
    }

    // Reverse the peptide index
    std::reverse(peptide_index.begin(), peptide_index.end());

    // Re-insert the fixed positions at their original places
    for (IndexType::iterator it = idx.begin(); it != idx.end(); ++it)
    {
      peptide_index.insert(peptide_index.begin() + *it, *it);
    }

    // Apply the reversed index to create the new sequence
    for (Size i = 0; i < peptide_index.size(); i++)
    {
      reversed[i] = sequence[peptide_index[i]];
    }

    // Relocate modifications based on index mapping
    for (auto& mod : reversed_mods)
    {
      // C and N terminal mods are implicitly not reversed (positions -1 and sequence.size())
      if (mod.location >= 0 && mod.location < static_cast<int>(sequence.size()))
      {
        for (Size k = 0; k < peptide_index.size(); k++)
        {
          if (static_cast<int>(peptide_index[k]) == mod.location)
          {
            mod.location = static_cast<int>(k);
            break;
          }
        }
      }
    }

    return std::make_pair(reversed, reversed_mods);
  }

  std::pair<std::string, std::vector<OpenSwath::LightModification>> MRMDecoy::pseudoreversePeptideLight_(
      const std::string& sequence,
      const std::vector<OpenSwath::LightModification>& modifications) const
  {
    // Match heavy version: pseudoreversePeptide_ calls reversePeptide with no const_pattern (defaults to empty)
    return MRMDecoy::reversePeptideLight(sequence, modifications, false, true);
  }

  void MRMDecoy::switchKRLight(std::string& sequence)
  {
    static std::string aa[] =
    {
      "A", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "M", "F", "S", "T", "W",
      "Y", "V"
    };
    constexpr int aa_size = 17;

    if (sequence.empty()) return;

    Size lastAA = sequence.size() - 1;
    if (sequence[lastAA] == 'K')
    {
      sequence[lastAA] = 'R';
    }
    else if (sequence[lastAA] == 'R')
    {
      sequence[lastAA] = 'K';
    }
    else
    {
      // Use portable FNV-1a hash of sequence to deterministically select replacement AA.
      // This ensures the same peptide always gets the same replacement regardless of
      // processing order, making results reproducible across platforms (including ARM64).
      uint64_t hash = 14695981039346656037ULL; // FNV offset basis
      for (char c : sequence)
      {
        hash ^= static_cast<uint64_t>(c);
        hash *= 1099511628211ULL; // FNV prime
      }
      int res_pos = static_cast<int>(hash % aa_size);
      sequence[lastAA] = aa[res_pos][0];
    }
  }

  bool MRMDecoy::hasCNterminalModsLight_(
      const std::vector<OpenSwath::LightModification>& modifications,
      size_t sequence_length,
      bool checkCterminalAA)
  {
    for (const auto& mod : modifications)
    {
      // N-terminal modification at position -1
      if (mod.location == -1)
      {
        return true;
      }
      // C-terminal modification at position sequence_length
      if (mod.location == static_cast<int>(sequence_length))
      {
        return true;
      }
      // Check C-terminal AA position if requested
      if (checkCterminalAA && mod.location == static_cast<int>(sequence_length) - 1)
      {
        return true;
      }
    }
    return false;
  }

  std::pair<std::string, std::vector<OpenSwath::LightModification>> MRMDecoy::shufflePeptideLight(
      const std::string& sequence,
      const std::vector<OpenSwath::LightModification>& modifications,
      const double identity_threshold,
      int seed,
      const int max_attempts) const
  {
    if (seed == -1)
    {
      seed = time(nullptr);
    }

    std::string shuffled = sequence;
    std::vector<OpenSwath::LightModification> shuffled_mods = modifications;

    // Working copy of input that may be mutated during iterations
    std::string peptide_seq = sequence;

    boost::mt19937 generator(seed);
    boost::uniform_int<> uni_dist;
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > pseudoRNG(generator, uni_dist);

    static std::string aa[] =
    {
      "A", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "M", "F", "S", "T", "W",
      "Y", "V"
    };
    constexpr int aa_size = 17;

    int attempts = 0;
    // Loop: attempt to shuffle and check whether difference is large enough
    // Note: compare against peptide_seq (not sequence) since peptide_seq may be mutated
    while (AASequenceIdentity(peptide_seq, shuffled) > identity_threshold &&
           attempts < max_attempts)
    {
      // Block tryptic residues and N-/C-terminus from shuffling
      MRMDecoy::IndexType idx = findFixedAndTermResidues_(peptide_seq);

      shuffled = peptide_seq;
      shuffled_mods = modifications;

      std::vector<Size> peptide_index;
      for (Size i = 0; i < peptide_seq.size(); i++)
      {
        peptide_index.push_back(i);
      }

      // Erase the indices where K/P/R are (from the back to preserve indices)
      for (IndexType::reverse_iterator it = idx.rbegin(); it != idx.rend(); ++it)
      {
        peptide_index.erase(peptide_index.begin() + *it);
      }

      // Shuffle the peptide index (Fisher-Yates shuffle)
      if (peptide_index.begin() != peptide_index.end())
      {
        for (std::vector<Size>::iterator pI_it = peptide_index.begin() + 1; pI_it != peptide_index.end(); ++pI_it)
        {
          std::iter_swap(pI_it, peptide_index.begin() + pseudoRNG((pI_it - peptide_index.begin()) + 1));
        }
      }

      // Re-insert the fixed positions at their original places
      for (IndexType::iterator it = idx.begin(); it != idx.end(); ++it)
      {
        peptide_index.insert(peptide_index.begin() + *it, *it);
      }

      // Apply shuffled index to create new sequence
      for (Size i = 0; i < peptide_index.size(); i++)
      {
        shuffled[i] = peptide_seq[peptide_index[i]];
      }

      // Relocate modifications based on index mapping
      for (Size j = 0; j < shuffled_mods.size(); j++)
      {
        // C and N terminal mods are implicitly not shuffled (positions -1 and sequence.size())
        if (shuffled_mods[j].location >= 0 && shuffled_mods[j].location < static_cast<int>(peptide_seq.size()))
        {
          for (Size k = 0; k < peptide_index.size(); k++)
          {
            if (static_cast<int>(peptide_index[k]) == shuffled_mods[j].location)
            {
              shuffled_mods[j].location = static_cast<int>(k);
              break;
            }
          }
        }
      }

      ++attempts;

      // Every 10 attempts (at 9, 19, 29...), mutate a random non-fixed AA to help convergence
      // This matches the heavy version's timing (attempts % 10 == 9)
      if (attempts % 10 == 9)
      {
        // Important: pick the new amino acid FIRST (matching heavy version's RNG call order)
        int res_pos = (pseudoRNG() % aa_size);

        // Find positions that are not modified and not N/C terminal
        int pep_pos = -1;
        size_t pos_trials = 0;
        while (pep_pos < 0 && pos_trials < shuffled.size())
        {
          pep_pos = (pseudoRNG() % shuffled.size());
          // Check if position is modified or is N/C terminus
          bool is_modified = false;
          for (const auto& mod : shuffled_mods)
          {
            if (mod.location == pep_pos)
            {
              is_modified = true;
              break;
            }
          }
          if (is_modified || (pep_pos == 0) || (pep_pos == static_cast<int>(shuffled.size() - 1)))
          {
            pep_pos = -1;
          }
          else
          {
            // Mutate the amino acid
            shuffled[pep_pos] = aa[res_pos][0];
          }
          ++pos_trials;
        }
        // Important: persist the mutation for next iteration (like heavy version)
        peptide_seq = shuffled;
      }
    }

    return std::make_pair(shuffled, shuffled_mods);
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
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " due to C/N-terminal modifications" << '\n';
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
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " due to C/N-terminal modifications" << '\n';
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
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " due to C/N-terminal modifications" << '\n';
          exclusion_peptides.insert(peptide.id);
        }
        else if (do_switchKR) { switchKR(peptide); }
      }

      // Check that the decoy precursor does not happen to be a target precursor AND is not already present.
      // Use getModifiedPeptideSequence_ to match the key format used when populating allPeptideSequences above.
      const std::string peptide_key = MRMDecoy::getModifiedPeptideSequence_(peptide) + String(peptide.getChargeState());
      if (allPeptideSequences.find(peptide_key) != allPeptideSequences.end())
      {
        OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " since decoy peptide is also a target peptide or this decoy peptide is already present" << '\n';
        exclusion_peptides.insert(peptide.id);
      }
      else
      {
        // Since this decoy will be added, add it to the precursor map so that the same decoy is not added twice
        OPENMS_LOG_DEBUG << "[peptide] adding " << peptide.id << " to master list of peptides " << '\n';
        allPeptideSequences[peptide_key] = peptide.id;
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
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << decoy_tr.getPeptideRef() << " due to missing annotation" << '\n';
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
        OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << " due to missing transitions" << '\n';
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
        OPENMS_LOG_DEBUG << "[protein] Skipping " << protein.id << " due to missing peptides" << '\n';
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
        OPENMS_LOG_DEBUG << "[peptide] Skipping " << decoy_compound.id << " - cannot parse sequence" << '\n';
        exclusion_peptides.insert(decoy_compound.id);
        continue;
      }

      // Use pure light-path decoy generation methods (no heavy Peptide allocation)
      std::string decoy_sequence;
      std::vector<OpenSwath::LightModification> decoy_mods;

      // Get unmodified sequence for decoy generation
      std::string unmodified_sequence = original_sequence.toUnmodifiedString();

      // Apply decoy method using light methods
      if (method == "pseudo-reverse")
      {
        if (hasCNterminalModsLight_(decoy_compound.modifications, unmodified_sequence.size(), do_switchKR))
        {
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << decoy_compound.id << " due to C/N-terminal modifications" << '\n';
          exclusion_peptides.insert(decoy_compound.id);
          continue;
        }
        auto result = pseudoreversePeptideLight_(unmodified_sequence, decoy_compound.modifications);
        decoy_sequence = result.first;
        decoy_mods = result.second;
        if (do_switchKR) { switchKRLight(decoy_sequence); }
      }
      else if (method == "reverse")
      {
        if (hasCNterminalModsLight_(decoy_compound.modifications, unmodified_sequence.size(), false))
        {
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << decoy_compound.id << " due to C/N-terminal modifications" << '\n';
          exclusion_peptides.insert(decoy_compound.id);
          continue;
        }
        // Note: heavy reversePeptide_() calls reversePeptide(peptide, false, false) with no const_pattern
        // so we must also use empty const_pattern here (the default)
        auto result = reversePeptideLight(unmodified_sequence, decoy_compound.modifications, false, false);
        decoy_sequence = result.first;
        decoy_mods = result.second;
      }
      else if (method == "shuffle")
      {
        auto result = shufflePeptideLight(unmodified_sequence, decoy_compound.modifications, identity_threshold, -1, max_attempts);
        decoy_sequence = result.first;
        decoy_mods = result.second;
        if (do_switchKR && hasCNterminalModsLight_(decoy_mods, decoy_sequence.size(), do_switchKR))
        {
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << decoy_compound.id << " due to C/N-terminal modifications" << '\n';
          exclusion_peptides.insert(decoy_compound.id);
          continue;
        }
        else if (do_switchKR)
        {
          switchKRLight(decoy_sequence);
        }
      }
      else
      {
        decoy_sequence = unmodified_sequence;
        decoy_mods = decoy_compound.modifications;
      }

      // Build modified sequence string with UniMod notation (must be done before duplicate check)
      // This matches the Heavy path which uses getModifiedPeptideSequence_() for duplicate detection
      String full_decoy_sequence;
      for (int loc = -1; loc <= static_cast<int>(decoy_sequence.size()); loc++)
      {
        if (loc > -1 && loc < static_cast<int>(decoy_sequence.size()))
        {
          full_decoy_sequence += decoy_sequence[loc];
        }
        // Add modifications at this location
        for (const auto& mod : decoy_mods)
        {
          if (mod.location == loc)
          {
            full_decoy_sequence += "(UniMod:" + String(mod.unimod_id) + ")";
          }
        }
      }
      decoy_compound.sequence = full_decoy_sequence;

      // Check for duplicates using MODIFIED sequence + charge (matching Heavy path behavior)
      // The Heavy path uses getModifiedPeptideSequence_(peptide) + charge for storing decoys,
      // which includes UniMod annotations. This ensures peptides with the same unmodified
      // sequence but different modifications are NOT considered duplicates.
      std::string decoy_key = full_decoy_sequence + String(decoy_compound.charge);
      if (allPeptideSequences.find(decoy_key) != allPeptideSequences.end())
      {
        OPENMS_LOG_DEBUG << "[peptide] Skipping " << decoy_compound.id << " since decoy peptide is also a target peptide or this decoy peptide is already present" << '\n';
        exclusion_peptides.insert(decoy_compound.id);
        continue;
      }
      // Add to map with modified sequence (matching Heavy path line 612)
      OPENMS_LOG_DEBUG << "[peptide] adding " << decoy_compound.id << " to master list of peptides" << '\n';
      allPeptideSequences[decoy_key] = decoy_compound.id;

      // Update modifications
      decoy_compound.modifications = decoy_mods;

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
          decoy_tr.setDecoy(true);
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
        OPENMS_LOG_DEBUG << "[transition] Skipping transitions for " << peptide_ref << " - cannot parse sequence" << '\n';
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
        if (!tr->isDetectingTransition() || tr->getDecoy())
        {
          continue;
        }

        OpenSwath::LightTransition decoy_tr = *tr;
        decoy_tr.transition_name = decoy_tag + tr->transition_name;
        decoy_tr.peptide_ref = decoy_peptide_ref;
        decoy_tr.setDecoy(true);
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
          // Update fragment type/number/charge from annotation (format: b4^1, y10^2, etc.)
          if (!targetion.first.empty())
          {
            decoy_tr.setFragmentType(targetion.first.substr(0, 1));
            std::string num_str;
            size_t i = 1;
            for (; i < targetion.first.size() && std::isdigit(targetion.first[i]); ++i)
            {
              num_str += targetion.first[i];
            }
            if (!num_str.empty())
            {
              decoy_tr.fragment_nr = static_cast<int16_t>(std::stoi(num_str));
            }
            // Extract charge after '^'
            if (i < targetion.first.size() && targetion.first[i] == '^')
            {
              std::string charge_str;
              for (size_t j = i + 1; j < targetion.first.size() && std::isdigit(targetion.first[j]); ++j)
              {
                charge_str += targetion.first[j];
              }
              if (!charge_str.empty())
              {
                decoy_tr.fragment_charge = static_cast<int8_t>(std::stoi(charge_str));
              }
            }
          }
          decoy_transitions.push_back(decoy_tr);
        }
        else
        {
          // Transition could not be annotated, exclude whole peptide (matching Heavy path behavior)
          OPENMS_LOG_DEBUG << "[peptide] Skipping " << decoy_peptide_ref << " due to missing annotation" << '\n';
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

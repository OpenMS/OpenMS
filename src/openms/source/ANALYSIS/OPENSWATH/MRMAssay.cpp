// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <boost/lexical_cast.hpp>
#include <regex>
#include <unordered_set>
#include <map>

using namespace std;

namespace OpenMS
{
  MRMAssay::MRMAssay() = default;

  MRMAssay::~MRMAssay() = default;

  std::vector<std::string> MRMAssay::getMatchingPeptidoforms_(const double fragment_ion,
                                                              const FragmentSeqMap& ions,
                                                              const double mz_threshold)
  {
    std::vector<std::string> isoforms;

    for (const auto& i_it : ions) // map: "fragment m/z" -> "modified sequence"
    {
      if (i_it.first - mz_threshold <= fragment_ion && i_it.first + mz_threshold >= fragment_ion)
      {
        isoforms.push_back(i_it.second);
      }
    }

    std::sort(isoforms.begin(), isoforms.end());
    isoforms.erase(std::unique(isoforms.begin(), isoforms.end()), isoforms.end());

    return isoforms;
  }

  int MRMAssay::getSwath_(const std::vector<std::pair<double, double> >& swathes, const double precursor_mz)
  {
    int swath = -1;

    // we go through all swaths in ascending order and if the transitions falls
    // in overlap, only the upper swath will be used and checked.
    for (auto it = swathes.begin(); it != swathes.end(); ++it)
    {
      if (precursor_mz >= it->first && precursor_mz <= it->second)
      {
        swath = it - swathes.begin();
      }
    }

    if (swath != -1)
    {
      return swath;
    }
    else
    {
      return -1;
    }
  }

  bool MRMAssay::isInSwath_(const std::vector<std::pair<double, double> >& swathes, const double precursor_mz, const double product_mz)
  {
    int swath_idx = getSwath_(swathes, precursor_mz);

    if (swath_idx == -1) { return true; } // remove all transitions that are not in swath range
    else
    {
      std::pair<double, double> swath = swathes[getSwath_(swathes, precursor_mz)];

      if (product_mz >= swath.first && product_mz <= swath.second)
      {
        return true;
      }
      else
      { 
        return false;
      }
    }
  }

  std::string MRMAssay::getRandomSequence_(size_t sequence_size, boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
                                           pseudoRNG)
  {
    std::string aa[] =
    {
      "A", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "M", "F", "S", "T", "W",
      "Y", "V"
    };
    size_t aa_size = 17;

    std::string peptide_sequence = "";

    for (size_t i = 0; i < sequence_size; ++i)
    {
      size_t pos = (pseudoRNG() % aa_size);
      peptide_sequence += aa[pos];
    }

    return peptide_sequence;
  }

  std::vector<std::vector<size_t> > MRMAssay::nchoosekcombinations_(const std::vector<size_t>& n, size_t k)
  {
    std::vector<std::vector<size_t> > combinations;

    std::string bitmask(k, 1);
    bitmask.resize(n.size(), 0);

    do
    {
      std::vector<size_t> combination;
      for (size_t i = 0; i < n.size(); ++i)
      {
        if (bitmask[i])
        {
          combination.push_back(n[i]);
        }
      }
      combinations.push_back(combination);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    return combinations;
  }

  std::vector<OpenMS::AASequence> MRMAssay::addModificationsSequences_(const std::vector<OpenMS::AASequence>& sequences, const std::vector<std::vector<size_t> >& mods_combs, const OpenMS::String& modification)
  {
    std::vector<OpenMS::AASequence> modified_sequences;
    bool multi_mod_switch = false;
    bool skip_invalid_mod_seq = false;

    OpenMS::ModificationsDB* ptr = ModificationsDB::getInstance();
    std::set<const ResidueModification*> modifiable_nterm;
    ptr->searchModifications(modifiable_nterm, modification, "", ResidueModification::N_TERM);
    std::set<const ResidueModification*> modifiable_cterm;
    ptr->searchModifications(modifiable_cterm, modification, "", ResidueModification::C_TERM);
    for (std::vector<OpenMS::AASequence>::const_iterator sq_it = sequences.begin(); sq_it != sequences.end(); ++sq_it)
    {
      for (std::vector<std::vector<size_t> >::const_iterator mc_it = mods_combs.begin(); mc_it != mods_combs.end(); ++mc_it)
      {
        multi_mod_switch = false;
        skip_invalid_mod_seq = false;
        OpenMS::AASequence temp_sequence = *sq_it;
        for (std::vector<size_t>::const_iterator pos_it = mc_it->begin(); pos_it != mc_it->end(); ++pos_it)
        {
          if (*pos_it == 0)
          {
            // Check first to make sure ending residue is NTerm modifiable
            if ( !modifiable_nterm.empty() && (temp_sequence[0].getOneLetterCode() == OpenMS::String((*modifiable_nterm.begin())->getOrigin()) || (*modifiable_nterm.begin())->getOrigin() == 'X') ) 
            {
              temp_sequence.setNTerminalModification(modification);
            } 
            else 
            {
              OPENMS_LOG_DEBUG << "[addModificationsSequences_] Skipping addition of N-Term " << OpenMS::String((*modifiable_nterm.begin())->getId()) <<
                                   " to last residue (" << temp_sequence[temp_sequence.size() - 1].getOneLetterCode() << ") of peptide " << temp_sequence.toUniModString() << 
                                   " , because it does not match viable N-Term residue specificity (" <<
                                   OpenMS::String((*modifiable_nterm.begin())->getOrigin()) << ") in ModificationDB." << std::endl;
              skip_invalid_mod_seq = true;
            }
          }
          else if (*pos_it == temp_sequence.size() + 1)
          {
            // Check first to make sure ending residue is CTerm modifiable
            if ( !modifiable_cterm.empty() && (temp_sequence.toUnmodifiedString().back() == (*modifiable_cterm.begin())->getOrigin() || (*modifiable_cterm.begin())->getOrigin() == 'X') )
            {
              temp_sequence.setCTerminalModification(modification);
            } 
            else 
            {
              OPENMS_LOG_DEBUG << "[addModificationsSequences_] Skipping addition of C-Term " << OpenMS::String((*modifiable_cterm.begin())->getId()) <<
                                   " to last residue (" << temp_sequence.toUnmodifiedString().back() << ") of peptide " << temp_sequence.toUniModString() << 
                                   " , because it does not match viable C-Term residue specificity (" <<
                                   OpenMS::String((*modifiable_cterm.begin())->getOrigin()) << ") in ModificationDB." << std::endl;
              skip_invalid_mod_seq = true;
            }
          }
          else
          {
            if (!temp_sequence[*pos_it - 1].isModified())
            {
              temp_sequence.setModification(*pos_it - 1, modification);
            }
            else
            {
              multi_mod_switch = true;
            }
          }
        }
        if (skip_invalid_mod_seq) { continue; }
        if (!multi_mod_switch) { modified_sequences.push_back(temp_sequence); }
      }
    }

    return modified_sequences;
  }

  std::vector<OpenMS::AASequence> MRMAssay::generateTheoreticalPeptidoforms_(const OpenMS::AASequence& sequence)
  {
    std::map<OpenMS::String, size_t> mods;
    std::vector<OpenMS::AASequence> sequences = {AASequence::fromString(sequence.toUnmodifiedString())};

    OpenMS::ModificationsDB* ptr = ModificationsDB::getInstance();

    if (sequence.hasNTerminalModification())
    {
      mods[sequence.getNTerminalModificationName()] += 1;
    }

    if (sequence.hasCTerminalModification())
    {
      mods[sequence.getCTerminalModificationName()] += 1;
    }

    for (size_t i = 0; i < sequence.size(); ++i)
    {
      if (sequence[i].isModified())
      {
        mods[sequence.getResidue(i).getModificationName()] += 1;
      }
    }

    // For each modification, create all (n choose k) theoretical peptidoforms
    for (const auto& mod_it : mods)
    {
      std::vector<size_t> mods_res;

      std::set<const ResidueModification*> modifiable_nterm;
      ptr->searchModifications(modifiable_nterm, mod_it.first, "", ResidueModification::N_TERM);
      if (!modifiable_nterm.empty())
      {
        mods_res.push_back(0);
      }

      std::set<const ResidueModification*> modifiable_cterm;
      ptr->searchModifications(modifiable_cterm, mod_it.first, "", ResidueModification::C_TERM);
      if (!modifiable_cterm.empty())
      {
        mods_res.push_back(sequence.size() + 1);
      }

      for (size_t i = 0; i < sequence.size(); ++i)
      {
        std::set<const ResidueModification*> modifiable_residues;
        ptr->searchModifications(modifiable_residues, mod_it.first, sequence.getResidue(i).getOneLetterCode(), ResidueModification::ANYWHERE);
        if (!modifiable_residues.empty())
        {
          mods_res.push_back(i + 1);
        }
      }
      std::vector<std::vector<size_t> > mods_combs = nchoosekcombinations_(mods_res, mod_it.second);
      sequences = addModificationsSequences_(sequences, mods_combs, mod_it.first);
    }
    return sequences;
  }

  std::vector<OpenMS::AASequence> MRMAssay::generateTheoreticalPeptidoformsDecoy_(const OpenMS::AASequence& sequence, const OpenMS::AASequence& decoy_sequence)
  {
    std::map<OpenMS::String, size_t> mods;
    std::vector<OpenMS::AASequence> decoy_sequences;
    decoy_sequences.push_back(AASequence::fromString(decoy_sequence.toUnmodifiedString()));

    OpenMS::ModificationsDB* ptr = ModificationsDB::getInstance();

    if (sequence.hasNTerminalModification())
    {
      mods[sequence.getNTerminalModificationName()] += 1;
    }

    if (sequence.hasCTerminalModification())
    {
      mods[sequence.getCTerminalModificationName()] += 1;
    }

    for (size_t i = 0; i < sequence.size(); ++i)
    {
      if (sequence[i].isModified())
      {
        mods[sequence.getResidue(i).getModificationName()] += 1;
      }
    }

    for (const auto& mod_it : mods)
    {
      std::vector<size_t> mods_res;

      std::set<const ResidueModification*> modifiable_nterm;
      ptr->searchModifications(modifiable_nterm, mod_it.first, "", ResidueModification::N_TERM);
      if (!modifiable_nterm.empty())
      {
        mods_res.push_back(0);
      }

      std::set<const ResidueModification*> modifiable_cterm;
      ptr->searchModifications(modifiable_cterm, mod_it.first, "", ResidueModification::C_TERM);
      if (!modifiable_cterm.empty())
      {
        mods_res.push_back(sequence.size() + 1);
      }

      for (size_t i = 0; i < sequence.size(); ++i)
      {
        std::set<const ResidueModification*> modifiable_residues;
        ptr->searchModifications(modifiable_residues, mod_it.first, sequence.getResidue(i).getOneLetterCode(), ResidueModification::ANYWHERE);
        if (!modifiable_residues.empty())
        {
          mods_res.push_back(i + 1);
        }
      }
      std::vector<std::vector<size_t> > mods_combs = nchoosekcombinations_(mods_res, mod_it.second);
      decoy_sequences = addModificationsSequences_(decoy_sequences, mods_combs, mod_it.first);
    }
    return decoy_sequences;
  }

  void MRMAssay::generateTargetInSilicoMap_(const OpenMS::TargetedExperiment& exp,
                                            const std::vector<String>& fragment_types,
                                            const std::vector<size_t>& fragment_charges,
                                            bool enable_specific_losses,
                                            bool enable_unspecific_losses,
                                            bool enable_ms2_precursors,
                                            const std::vector<std::pair<double, double> >& swathes,
                                            int round_decPow,
                                            size_t max_num_alternative_localizations,
                                            SequenceMapT & TargetSequenceMap,
                                            IonMapT & TargetIonMap,
                                            PeptideMapT& TargetPeptideMap)
  {
    OpenMS::MRMIonSeries mrmis;

    // Step 1: Generate target in silico peptide map containing theoretical transitions
    Size progress = 0;
    startProgress(0, exp.getPeptides().size(), "Generation of target in silico peptide map");
    for (size_t i = 0; i < exp.getPeptides().size(); ++i)
    {
      setProgress(progress++);

      TargetedExperiment::Peptide peptide = exp.getPeptides()[i];
      OpenMS::AASequence peptide_sequence = TargetedExperimentHelper::getAASequence(peptide);
      int precursor_charge = 1;
      if (peptide.hasCharge()) 
      {
        precursor_charge = peptide.getChargeState();
      }
      double precursor_mz = peptide_sequence.getMZ(precursor_charge);
      int precursor_swath = getSwath_(swathes, precursor_mz);

      // Compute all alternative peptidoforms compatible with ModificationsDB
      const vector<AASequence> alternative_peptide_sequences = generateTheoreticalPeptidoforms_(peptide_sequence);  

      // Some permutations might be too complex, skip if threshold is reached
      if (alternative_peptide_sequences.size() > max_num_alternative_localizations)
      {
        OPENMS_LOG_DEBUG << "[uis] Peptide skipped (too many permutations possible): " << peptide.id << std::endl;
        continue;
      }

      // Iterate over all peptidoforms
      for (const auto& alt_aa : alternative_peptide_sequences)
      { 
        // Append peptidoform to index
        TargetSequenceMap[precursor_swath][alt_aa.toUnmodifiedString()].insert(alt_aa.toString());
        // Generate theoretical ion series
        auto ionseries = mrmis.getIonSeries(alt_aa, precursor_charge,
            fragment_types, fragment_charges, enable_specific_losses,
            enable_unspecific_losses);

        if (enable_ms2_precursors)
        {
          // Add precursor to theoretical transitions
          double prec_mz = Math::roundDecimal(precursor_mz, round_decPow);
          TargetIonMap[precursor_swath][alt_aa.toUnmodifiedString()].emplace_back(prec_mz, alt_aa.toString());
          TargetPeptideMap[peptide.id].emplace_back("MS2_Precursor_i0", prec_mz);
        }

        // Iterate over all theoretical transitions
        for (const auto& im_it : ionseries)
        {
          // Append transition to indices to find interfering transitions
          double fragment_mz = Math::roundDecimal(im_it.second, round_decPow);
          TargetIonMap[precursor_swath][alt_aa.toUnmodifiedString()].emplace_back(fragment_mz, alt_aa.toString());
          TargetPeptideMap[peptide.id].emplace_back(im_it.first, fragment_mz);
        }
      }
    }
    endProgress();
  }

  void MRMAssay::generateDecoySequences_(const SequenceMapT& TargetSequenceMap,
                                         std::map<String, String>& DecoySequenceMap, int shuffle_seed)
  {
    // Step 2a: Generate decoy sequences that share peptidoform properties with targets
    if (shuffle_seed == -1)
    {
      shuffle_seed = time(nullptr);
    }

    boost::mt19937 generator(shuffle_seed);
    boost::uniform_int<> uni_dist;
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > pseudoRNG(generator, uni_dist);

    Size progress = 0;
    startProgress(0, TargetSequenceMap.size(), "Target-decoy mapping");
    std::string decoy_peptide_string;

    // Iterate over swathes
    for (const auto& sm_it : TargetSequenceMap)
    {
      setProgress(progress++);
      // Iterate over each unmodified peptide sequence in current SWATH
      for (const auto& ta_it : sm_it.second)
      {
        // Get a random unmodified peptide sequence as base for later modification
        if (DecoySequenceMap[ta_it.first].empty())
        {
          decoy_peptide_string = getRandomSequence_(ta_it.first.size(), pseudoRNG);
        }
        else
        {
          decoy_peptide_string = DecoySequenceMap[ta_it.first];
        }

        // Iterate over all target peptidoforms for the current unmodified
        // peptide sequence and replace decoy residues with modified target
        // residues
        for (const auto & se_it : ta_it.second)
        {
          OpenMS::AASequence seq = AASequence::fromString(se_it);

          if (seq.hasNTerminalModification())
          {
            decoy_peptide_string = decoy_peptide_string.replace(0, 1, seq.getSubsequence(0, 1).toUnmodifiedString());
          }

          if (seq.hasCTerminalModification())
          {
            decoy_peptide_string = decoy_peptide_string.replace(decoy_peptide_string.size() - 1, 1, seq.getSubsequence(decoy_peptide_string.size() - 1, 1).toUnmodifiedString());
          }

          for (size_t i = 0; i < seq.size(); ++i)
          {
            if (seq[i].isModified())
            {
              decoy_peptide_string = decoy_peptide_string.replace(i, 1, seq.getSubsequence(i, 1).toUnmodifiedString());
            }
          }
          DecoySequenceMap[ta_it.first] = decoy_peptide_string;
        }
      }
    }
    endProgress();
  }

  void MRMAssay::generateDecoyInSilicoMap_(const OpenMS::TargetedExperiment& exp,
                                           const std::vector<String>& fragment_types,
                                           const std::vector<size_t>& fragment_charges,
                                           bool enable_specific_losses,
                                           bool enable_unspecific_losses,
                                           bool enable_ms2_precursors,
                                           const std::vector<std::pair<double, double> >& swathes,
                                           int round_decPow,
                                           TargetDecoyMapT& TargetDecoyMap,
                                           PeptideMapT& TargetPeptideMap,
                                           std::map<String, String>& DecoySequenceMap,
                                           IonMapT & DecoyIonMap,
                                           PeptideMapT& DecoyPeptideMap)
  {
    MRMIonSeries mrmis;

    // Step 2b: Generate decoy in silico peptide map containing theoretical transitions
    Size progress = 0;
    startProgress(0, exp.getPeptides().size(), "Generation of decoy in silico peptide map");
    for (size_t i = 0; i < exp.getPeptides().size(); ++i)
    {
      setProgress(progress++);

      TargetedExperiment::Peptide peptide = exp.getPeptides()[i];
      int precursor_charge = 1;
      if (peptide.hasCharge()) 
      {
        precursor_charge = peptide.getChargeState();
      }

      // Skip if target peptide is not in map, e.g. permutation threshold was reached
      if (TargetPeptideMap.find(peptide.id) == TargetPeptideMap.end())
      {
        continue;
      }

      OpenMS::AASequence peptide_sequence = TargetedExperimentHelper::getAASequence(peptide);
      double precursor_mz = peptide_sequence.getMZ(precursor_charge);
      int precursor_swath = getSwath_(swathes, precursor_mz);

      // Copy properties of target peptide to decoy and get sequence from map
      TargetedExperiment::Peptide decoy_peptide = peptide;
      decoy_peptide.sequence = DecoySequenceMap[peptide.sequence];

      TargetDecoyMap[peptide.id] = decoy_peptide;
      OpenMS::AASequence decoy_peptide_sequence = TargetedExperimentHelper::getAASequence(decoy_peptide);

      // Compute all alternative peptidoforms compatible with ModificationsDB
      // Infers residue specificity from target sequence but is applied to decoy sequence
      const vector<AASequence> alternative_decoy_peptide_sequences = generateTheoreticalPeptidoformsDecoy_(peptide_sequence, decoy_peptide_sequence);

      // Iterate over all peptidoforms
      for (const auto& alt_aa : alternative_decoy_peptide_sequences)
      {
        // Generate theoretical ion series
        MRMIonSeries::IonSeries ionseries = mrmis.getIonSeries(alt_aa, precursor_charge, // use same charge state as target
                                                               fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses);  

        if (enable_ms2_precursors)
        {
          // Add precursor to theoretical transitions
          double prec_mz = Math::roundDecimal(precursor_mz, round_decPow);
          DecoyIonMap[precursor_swath][alt_aa.toUnmodifiedString()].emplace_back(prec_mz, alt_aa.toString());
          DecoyPeptideMap[peptide.id].emplace_back("MS2_Precursor_i0", prec_mz);
        }

        // Iterate over all theoretical transitions
        for (const auto& im_it : ionseries)
        {
          // Append transition to indices to find interfering transitions
          double fragment_mz = Math::roundDecimal(im_it.second, round_decPow);
          DecoyIonMap[precursor_swath][alt_aa.toUnmodifiedString()].emplace_back(fragment_mz, alt_aa.toString());
          DecoyPeptideMap[decoy_peptide.id].emplace_back(im_it.first, fragment_mz);
        }
      }
    }
    endProgress();
  }

 void MRMAssay::generateTargetAssays_(const OpenMS::TargetedExperiment& exp,
                                      TransitionVectorType& transitions,
                                      double mz_threshold,
                                      const std::vector<std::pair<double, double> >& swathes,
                                      int round_decPow,
                                      const PeptideMapT& TargetPeptideMap,
                                      const IonMapT & TargetIonMap)
  {
    MRMIonSeries mrmis;

    // Step 3: Generate target identification transitions
    Size progress = 0;
    startProgress(0, TargetPeptideMap.size(), "Generation of target identification transitions");

    // Iterate over all target peptides
    int transition_index = 0;
    for (const auto& pep_it : TargetPeptideMap)
    { 
      setProgress(progress++);

      TargetedExperiment::Peptide peptide = exp.getPeptideByRef(pep_it.first);
      int precursor_charge = 1;
      if (peptide.hasCharge()) 
      {
        precursor_charge = peptide.getChargeState();
      }
      AASequence peptide_sequence = TargetedExperimentHelper::getAASequence(peptide);
      int target_precursor_swath = getSwath_(swathes, peptide_sequence.getMZ(precursor_charge));

      // Sort all transitions and make them unique
      auto transition_vector = pep_it.second;
      std::sort(transition_vector.begin(), transition_vector.end());
      auto tr_vec_end = std::unique(transition_vector.begin(), transition_vector.end());

      // Iterate over all transitions
      for (auto tr_it = transition_vector.begin(); tr_it != tr_vec_end; ++tr_it)
      { 
        // Compute the set of peptidoforms mapping to this transition
        vector<string> isoforms = getMatchingPeptidoforms_(
            tr_it->second, TargetIonMap.at(target_precursor_swath).at(peptide_sequence.toUnmodifiedString()), mz_threshold);

        // Check that transition maps to at least one peptidoform
        if (!isoforms.empty())
        {
          ReactionMonitoringTransition trn;
          trn.setDetectingTransition(false);
          trn.setMetaValue("insilico_transition", "true");
          trn.setPrecursorMZ(Math::roundDecimal(peptide_sequence.getMZ(precursor_charge), round_decPow));
          trn.setProductMZ(tr_it->second);
          trn.setPeptideRef(peptide.id);
          mrmis.annotateTransitionCV(trn, tr_it->first);
          trn.setIdentifyingTransition(true);
          trn.setQuantifyingTransition(false);

          // Set transition name containing mapping to peptidoforms with potential peptidoforms enumerated in brackets
          String identifier = String(transition_index) + "_" + String("UIS") +  \
            "_{" + ListUtils::concatenate(isoforms, "|") + "}_" +  \
            String(trn.getPrecursorMZ()) + "_" + String(trn.getProductMZ()) + "_" + 
            String(peptide.getRetentionTime()) + "_" + tr_it->first; 
          trn.setName(identifier); 
          trn.setNativeID(identifier);
          trn.setMetaValue("Peptidoforms", ListUtils::concatenate(isoforms, "|"));

          OPENMS_LOG_DEBUG << "[uis] Transition " << trn.getNativeID() << std::endl;

          // Append transition
          transitions.push_back(trn);
        }
        transition_index++;
      }
      OPENMS_LOG_DEBUG << "[uis] Peptide " << peptide.id << std::endl;
    }
    endProgress();
  }

 void MRMAssay::generateDecoyAssays_(const OpenMS::TargetedExperiment& exp,
                                     TransitionVectorType& transitions,
                                     double mz_threshold,
                                     const std::vector<std::pair<double, double> >& swathes,
                                     int round_decPow,
                                     const PeptideMapT& DecoyPeptideMap,
                                     std::map<String, TargetedExperiment::Peptide>& TargetDecoyMap,
                                     const IonMapT& DecoyIonMap,
                                     const IonMapT& TargetIonMap)
  {
    MRMIonSeries mrmis;

    // Step 4: Generate decoy identification transitions
    Size progress = 0;
    startProgress(0, DecoyPeptideMap.size(), "Generation of decoy identification transitions");

    // Iterate over all decoy peptides
    int transition_index = 0;
    for (const auto & decoy_pep_it : DecoyPeptideMap)
    {
      setProgress(progress++);
      const TargetedExperiment::Peptide& target_peptide = exp.getPeptideByRef(decoy_pep_it.first);
      int precursor_charge = 1;
      if (target_peptide.hasCharge()) 
      {
        precursor_charge = target_peptide.getChargeState();
      }
      AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);
      int target_precursor_swath = getSwath_(swathes, target_peptide_sequence.getMZ(precursor_charge));

      TargetedExperiment::Peptide decoy_peptide = TargetDecoyMap[decoy_pep_it.first];
      OpenMS::AASequence decoy_peptide_sequence = TargetedExperimentHelper::getAASequence(decoy_peptide);

      // Sort all transitions and make them unique
      auto transition_vector = decoy_pep_it.second;
      std::sort(transition_vector.begin(), transition_vector.end());
      auto tr_vec_end = std::unique(transition_vector.begin(), transition_vector.end());

      // Iterate over all transitions
      for (auto decoy_tr_it = transition_vector.begin(); decoy_tr_it != tr_vec_end; ++decoy_tr_it)
      {
        // Check mapping of transitions to other peptidoforms
        vector<string> decoy_isoforms = getMatchingPeptidoforms_(
            decoy_tr_it->second, DecoyIonMap.at(target_precursor_swath).at(decoy_peptide_sequence.toUnmodifiedString()), mz_threshold);

        // Check that transition maps to at least one peptidoform
        if (!decoy_isoforms.empty())
        {
          ReactionMonitoringTransition trn;
          trn.setDecoyTransitionType(ReactionMonitoringTransition::DECOY);
          trn.setDetectingTransition(false);
          trn.setMetaValue("insilico_transition", "true");
          trn.setPrecursorMZ(Math::roundDecimal(target_peptide_sequence.getMZ(precursor_charge), round_decPow));
          trn.setProductMZ(decoy_tr_it->second);
          trn.setPeptideRef(decoy_peptide.id);
          mrmis.annotateTransitionCV(trn, decoy_tr_it->first);
          trn.setIdentifyingTransition(true);
          trn.setQuantifyingTransition(false);

          // Set transition name containing mapping to peptidoforms with potential peptidoforms enumerated in brackets
          String identifier = String(transition_index) + "_" + String("UISDECOY") +
                "_{" + ListUtils::concatenate(decoy_isoforms, "|") + "}_" +
                String(trn.getPrecursorMZ()) + "_" + String(trn.getProductMZ()) + "_" +
                String(decoy_peptide.getRetentionTime()) + "_" + decoy_tr_it->first;
          trn.setName(identifier); 
          trn.setNativeID(identifier);
          trn.setMetaValue("Peptidoforms", ListUtils::concatenate(decoy_isoforms, "|"));

          OPENMS_LOG_DEBUG << "[uis] Decoy transition " << trn.getNativeID() << std::endl;

          // Check if decoy transition is overlapping with target transition
          vector<string> target_isoforms_overlap = getMatchingPeptidoforms_(
              decoy_tr_it->second, TargetIonMap.at(target_precursor_swath).at(target_peptide_sequence.toUnmodifiedString()), mz_threshold);

          if (!target_isoforms_overlap.empty())
          {
            OPENMS_LOG_DEBUG << "[uis] Skipping overlapping decoy transition " << trn.getNativeID() << std::endl;
            continue;
          }
          else
          {
            // Append transition
            transitions.push_back(trn);
          }
        }
        transition_index++;
      }
    }
    endProgress();
  }

  void MRMAssay::reannotateTransitions(OpenMS::TargetedExperiment& exp,
                             double precursor_mz_threshold,
                             double product_mz_threshold,
                             const std::vector<String>& fragment_types,
                             const std::vector<size_t>& fragment_charges,
                             bool enable_specific_losses,
                             bool enable_unspecific_losses,
                             int round_decPow)
  {
    PeptideVectorType peptides;
    ProteinVectorType proteins;
    TransitionVectorType transitions;

    OpenMS::MRMIonSeries mrmis;

    // hash of the peptide reference containing all transitions
    MRMAssay::PeptideTransitionMapType peptide_trans_map;
    for (Size i = 0; i < exp.getTransitions().size(); i++)
    {
      peptide_trans_map[exp.getTransitions()[i].getPeptideRef()].push_back(&exp.getTransitions()[i]);
    }

    Size progress = 0;
    startProgress(0, exp.getTransitions().size(), "Annotating transitions");
    for (MRMAssay::PeptideTransitionMapType::iterator pep_it = peptide_trans_map.begin();
         pep_it != peptide_trans_map.end(); ++pep_it)
    {
      String peptide_ref = pep_it->first;

      TargetedExperiment::Peptide target_peptide = exp.getPeptideByRef(peptide_ref);
      OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);

      int precursor_charge = 1;
      if (target_peptide.hasCharge()) {precursor_charge = target_peptide.getChargeState();}

      MRMIonSeries::IonSeries target_ionseries = mrmis.getIonSeries(
                                                    target_peptide_sequence, precursor_charge, fragment_types,
                                                    fragment_charges, enable_specific_losses,
                                                    enable_unspecific_losses, round_decPow);

      // Generate theoretical precursor m.z
      double precursor_mz = target_peptide_sequence.getMZ(precursor_charge);
      precursor_mz = Math::roundDecimal(precursor_mz, round_decPow);

      for (Size i = 0; i < pep_it->second.size(); i++)
      {
        setProgress(++progress);
        ReactionMonitoringTransition tr = *(pep_it->second[i]);

        // Annotate transition from theoretical ion series
        std::pair<String, double> targetion = mrmis.annotateIon(target_ionseries, tr.getProductMZ(), product_mz_threshold);

        // Ensure that precursor m/z is within threshold
        if (std::fabs(tr.getPrecursorMZ() - precursor_mz) > precursor_mz_threshold)
        {
          targetion.first = "unannotated";
        }

        // Set precursor m/z to theoretical value
        tr.setPrecursorMZ(precursor_mz);

        // Set product m/z to theoretical value
        tr.setProductMZ(targetion.second);

        // Skip unannotated transitions from previous step
        if (targetion.first == "unannotated")
        {
          OPENMS_LOG_DEBUG << "[unannotated] Skipping " << target_peptide_sequence.toString() 
            << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() 
            << " " << tr.getMetaValue("annotation") << std::endl;
          continue;
        }
        else
        {
          OPENMS_LOG_DEBUG << "[selected] " << target_peptide_sequence.toString() << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() << " " << tr.getMetaValue("annotation") << std::endl;
        }

        // Set CV terms
        mrmis.annotateTransitionCV(tr, targetion.first);

        // Add reference to parent precursor
        tr.setPeptideRef(target_peptide.id);

        // Append transition
        transitions.push_back(tr);
      }
    }
    endProgress();

    exp.setTransitions(std::move(transitions));
  }

  void MRMAssay::restrictTransitions(OpenMS::TargetedExperiment& exp, double lower_mz_limit, double upper_mz_limit, const std::vector<std::pair<double, double> >& swathes)
  {
    OpenMS::MRMIonSeries mrmis;
    PeptideVectorType peptides;
    ProteinVectorType proteins;
    TransitionVectorType transitions;

    Size progress = 0;
    startProgress(0, exp.getTransitions().size(), "Restricting transitions");
    for (Size i = 0; i < exp.getTransitions().size(); ++i)
    {
      setProgress(++progress);
      ReactionMonitoringTransition tr = exp.getTransitions()[i];

      const TargetedExperiment::Peptide& target_peptide = exp.getPeptideByRef(tr.getPeptideRef());
      OpenMS::AASequence target_peptide_sequence = TargetedExperimentHelper::getAASequence(target_peptide);

      // Check annotation for unannotated interpretations
      if (!tr.getProduct().getInterpretationList().empty())
      {
        // Check if transition is unannotated at primary annotation and if yes, skip
        if (tr.getProduct().getInterpretationList()[0].iontype == TargetedExperiment::IonType::NonIdentified)
        {
          OPENMS_LOG_DEBUG << "[unannotated] Skipping " << target_peptide_sequence 
            << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() 
            << " " << tr.getMetaValue("annotation") << std::endl;
          continue;
        }
      }

      // Check if product m/z falls into swath from precursor m/z and if yes, skip
      if (!swathes.empty())
      {
        if (MRMAssay::isInSwath_(swathes, tr.getPrecursorMZ(), tr.getProductMZ()))
        {
          OPENMS_LOG_DEBUG << "[swath] Skipping " << target_peptide_sequence << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() << std::endl;
          continue;
        }
      }

      // Check if product m/z is outside of m/z boundaries and if yes, skip
      if (tr.getProductMZ() < lower_mz_limit || tr.getProductMZ() > upper_mz_limit)
      {
        OPENMS_LOG_DEBUG << "[mz_limit] Skipping " << target_peptide_sequence << " PrecursorMZ: " << tr.getPrecursorMZ() << " ProductMZ: " << tr.getProductMZ() << std::endl;
        continue;
      }

      // Append transition
      transitions.push_back(tr);
    }

    exp.setTransitions(std::move(transitions));
    endProgress();
  }

  void MRMAssay::detectingTransitions(OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions)
  {
    PeptideVectorType peptides;
    ProteinVectorType proteins;
    TransitionVectorType transitions;

    std::unordered_set<String> peptide_ids;
    std::unordered_set<String> ProteinList;

    std::map<String, TransitionVectorType> TransitionsMap;

    // Generate a map of peptides to transitions for easy access
    for (Size i = 0; i < exp.getTransitions().size(); ++i)
    {
      ReactionMonitoringTransition tr = exp.getTransitions()[i];

      if (TransitionsMap.find(tr.getPeptideRef()) == TransitionsMap.end())
      {
        TransitionsMap[tr.getPeptideRef()];
      }

      TransitionsMap[tr.getPeptideRef()].push_back(tr);
    }

    Size progress = 0;
    startProgress(0, TransitionsMap.size() + exp.getPeptides().size() + exp.getProteins().size(), "Select detecting transitions");
    for (std::map<String, TransitionVectorType>::iterator m = TransitionsMap.begin();
         m != TransitionsMap.end(); ++m)
    {
      setProgress(++progress);
      // Ensure that all precursors have the minimum number of transitions
      if (m->second.size() >= (Size)min_transitions)
      {
        // LibraryIntensity stores all reference transition intensities of a precursor
        std::vector<double> LibraryIntensity;
        for (TransitionVectorType::iterator tr_it = m->second.begin(); tr_it != m->second.end(); ++tr_it)
        {
          LibraryIntensity.push_back(boost::lexical_cast<double>(tr_it->getLibraryIntensity()));
        }

        // Sort by intensity, reverse and delete all elements after max_transitions to find the best candidates
        std::sort(LibraryIntensity.begin(), LibraryIntensity.end());
        std::reverse(LibraryIntensity.begin(), LibraryIntensity.end());
        if ((Size)max_transitions < LibraryIntensity.size())
        {
          std::vector<double>::iterator start_delete = LibraryIntensity.begin();
          std::advance(start_delete, max_transitions);
          LibraryIntensity.erase(start_delete, LibraryIntensity.end());
        }

        // Check if transitions are among the ones with maximum intensity
        // If several transitions have the same intensities ensure restriction max_transitions
        Size j = 0; // transition number index
        for (TransitionVectorType::iterator tr_it = m->second.begin(); tr_it != m->second.end(); ++tr_it)
        {
          ReactionMonitoringTransition tr = *tr_it;

          if (
              (std::find(LibraryIntensity.begin(), LibraryIntensity.end(), boost::lexical_cast<double>(tr.getLibraryIntensity())) != LibraryIntensity.end()) &&
               tr.getDecoyTransitionType() != ReactionMonitoringTransition::DECOY &&
               j < (Size)max_transitions)
          {
            // Set meta value tag for detecting transition
            tr.setDetectingTransition(true);
            j += 1;
          }
          else
          {
            continue;
          }

          // Append transition
          transitions.push_back(tr);

          // Append transition_group_id to index
          peptide_ids.insert(tr.getPeptideRef());
        }
      }
    }

    for (const auto& peptide : exp.getPeptides())
    {
      setProgress(++progress);

      // Check if peptide has any transitions left
      if (peptide_ids.find(peptide.id) != peptide_ids.end())
      {
        peptides.push_back(peptide);
        for (const auto& protein_ref : peptide.protein_refs)
        {
          ProteinList.insert(protein_ref);
        }
      }
      else
      {
        OPENMS_LOG_DEBUG << "[peptide] Skipping " << peptide.id << std::endl;
      }
    }

    for (const auto& protein : exp.getProteins())
    {
      setProgress(++progress);

      // Check if protein has any peptides left
      if (ProteinList.find(protein.id) != ProteinList.end())
      {
        proteins.push_back(protein);
      }
      else
      {
        OPENMS_LOG_DEBUG << "[protein] Skipping " << protein.id << std::endl;
      }
    }

    exp.setTransitions(std::move(transitions));
    exp.setPeptides(std::move(peptides));
    exp.setProteins(std::move(proteins));

    endProgress();
  }

  void MRMAssay::uisTransitions(OpenMS::TargetedExperiment& exp,
                      const std::vector<String>& fragment_types,
                      const std::vector<size_t>& fragment_charges,
                      bool enable_specific_losses,
                      bool enable_unspecific_losses,
                      bool enable_ms2_precursors,
                      double mz_threshold,
                      const std::vector<std::pair<double, double> >& swathes,
                      int round_decPow,
                      size_t max_num_alternative_localizations,
                      int shuffle_seed,
                      bool disable_decoy_transitions)
  {
    OpenMS::MRMIonSeries mrmis;

    TransitionVectorType transitions = exp.getTransitions();

    // Different maps to store temporary data for fast access
    // TargetIonMap & DecoyIonMap: Store product m/z of all peptidoforms to find interfering transitions
    IonMapT TargetIonMap, DecoyIonMap;
    // TargetPeptideMap & DecoyPeptideMap: Store theoretical transitions of all peptidoforms
    PeptideMapT TargetPeptideMap, DecoyPeptideMap;
    // TargetSequenceMap, DecoySequenceMap & TargetDecoyMap: Link targets and UIS decoys
    SequenceMapT TargetSequenceMap;
    std::map<String, String> DecoySequenceMap;
    std::map<String, TargetedExperiment::Peptide> TargetDecoyMap;

    // Step 1: Generate target in silico peptide map containing theoretical transitions
    generateTargetInSilicoMap_(exp, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, enable_ms2_precursors, swathes, round_decPow, max_num_alternative_localizations, TargetSequenceMap, TargetIonMap, TargetPeptideMap);

    // Step 2: Generate target identification transitions
    generateTargetAssays_(exp, transitions, mz_threshold, swathes, round_decPow, TargetPeptideMap, TargetIonMap);

    if (!disable_decoy_transitions)
    {
      // Step 3a: Generate decoy sequences that share peptidoform properties with targets
      generateDecoySequences_(TargetSequenceMap, DecoySequenceMap, shuffle_seed);

      // Step 3b: Generate decoy in silico peptide map containing theoretical transitions
      generateDecoyInSilicoMap_(exp, fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses, enable_ms2_precursors, swathes, round_decPow, TargetDecoyMap, TargetPeptideMap, DecoySequenceMap, DecoyIonMap, DecoyPeptideMap);

      // Step 4: Generate decoy identification transitions
      generateDecoyAssays_(exp, transitions, mz_threshold, swathes, round_decPow, DecoyPeptideMap, TargetDecoyMap, DecoyIonMap, TargetIonMap);
    }

    exp.setTransitions(transitions);
  }

  void MRMAssay::filterMinMaxTransitionsCompound(OpenMS::TargetedExperiment& exp, int min_transitions, int max_transitions)
  {
    CompoundVectorType compounds;
    std::vector<String> compound_ids;
    TransitionVectorType transitions;

    std::map<String, TransitionVectorType> TransitionsMap;

    // Generate a map of compounds to transitions for easy access
    for (Size i = 0; i < exp.getTransitions().size(); ++i)
    {
      ReactionMonitoringTransition tr = exp.getTransitions()[i];

      if (TransitionsMap.find(tr.getCompoundRef()) == TransitionsMap.end())
      {
        TransitionsMap[tr.getCompoundRef()];
      }

      TransitionsMap[tr.getCompoundRef()].push_back(tr);
    }

    for (std::map<String, TransitionVectorType>::iterator m = TransitionsMap.begin();
         m != TransitionsMap.end(); ++m)
    {
        // Ensure that all precursors have the minimum number of transitions or are a decoy transitions
        if (m->second.size() >= (Size)min_transitions || m->second[0].getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
        {
        // LibraryIntensity stores all reference transition intensities of a precursor
        std::vector<double> LibraryIntensity;
        for (TransitionVectorType::iterator tr_it = m->second.begin(); tr_it != m->second.end(); ++tr_it)
        {
          LibraryIntensity.push_back(boost::lexical_cast<double>(tr_it->getLibraryIntensity()));
        }

        // Reverse-sort by intensity and delete all elements after max_transitions to find the best candidates
        std::sort(LibraryIntensity.begin(), LibraryIntensity.end(), std::greater<>());
        if ((Size)max_transitions < LibraryIntensity.size())
        {
          LibraryIntensity.resize(max_transitions);
        }

        // Check if transitions are among the ones with maximum intensity
        // If several transitions have the same intensities ensure restriction max_transitions
        Size j = 0; // transition number index
        for (TransitionVectorType::iterator tr_it = m->second.begin(); tr_it != m->second.end(); ++tr_it)
        {
          ReactionMonitoringTransition tr = *tr_it;
          if ((std::find(LibraryIntensity.begin(), LibraryIntensity.end(), boost::lexical_cast<double>(tr.getLibraryIntensity())) != LibraryIntensity.end()) && j < (Size)max_transitions)
          {
            // Set meta value tag for detecting transition
            tr.setDetectingTransition(true);
            j += 1;
          }
          else
          {
            continue;
          }

          // Append transition
          transitions.push_back(tr);

          // Append transition_group_id to index
          if (std::find(compound_ids.begin(), compound_ids.end(), tr.getCompoundRef()) == compound_ids.end())
          {
            compound_ids.push_back(tr.getCompoundRef());
          }
        }
      }
    }

    for (Size i = 0; i < exp.getCompounds().size(); ++i)
    {
      TargetedExperiment::Compound compound = exp.getCompounds()[i];

      // Check if compound has any transitions left
      if (std::find(compound_ids.begin(), compound_ids.end(), compound.id) != compound_ids.end())
      {
        compounds.push_back(compound);
      }
      else
      {
        OPENMS_LOG_DEBUG << "[compound] Skipping " << compound.id << " - not enough transitions."<< std::endl;
      }
    }
    exp.setTransitions(transitions);
    exp.setCompounds(compounds);
  }

  void MRMAssay::filterUnreferencedDecoysCompound(OpenMS::TargetedExperiment &exp)
  {
    vector<TargetedExperiment::Compound> compounds = exp.getCompounds();
    vector<std::string> descriptions_targets;
    vector<std::string> descriptions_decoys;
    vector<std::string> difference_target_decoys;
    vector<std::pair<std::string, std::string>> reference_decoys;
    vector<std::string> single_decoy_id;
    String decoy_suffix = "_decoy";

    for (const auto &it : compounds)
    {
      // extract potential target TransitionIds based on the decoy annotation '0_CompoundName_decoy_[M+H]+_448_0'
      if (it.id.find("decoy") != std::string::npos)
      {
        String current_decoy = it.id;
        String potential_target = current_decoy;
        potential_target.erase(potential_target.find(decoy_suffix), decoy_suffix.size());
        descriptions_decoys.emplace_back(potential_target);
        reference_decoys.emplace_back(std::make_pair(current_decoy, potential_target));
      }
      else
      {
        descriptions_targets.emplace_back(it.id);
      }
    }

    // compare the actual target TransitionsIds with the potential ones
    difference_target_decoys.clear();
    std::sort(descriptions_targets.begin(), descriptions_targets.end());
    std::sort(descriptions_decoys.begin(), descriptions_decoys.end());
    std::set_difference(descriptions_decoys.begin(),
                        descriptions_decoys.end(),
                        descriptions_targets.begin(),
                        descriptions_targets.end(),
                        std::inserter(difference_target_decoys, difference_target_decoys.begin()));

    // translate the potential targets back to decoy annotations
    for (const auto &it : difference_target_decoys)
    {
      auto iter = std::find_if(reference_decoys.begin(),
                               reference_decoys.end(),
                               [&it](const std::pair<String, String> &element) { return element.second == it; });
      if (iter != reference_decoys.end())
      {
        single_decoy_id.emplace_back(iter->first);
      }
    }

    // remove decoy compound due to missing target
    vector<TargetedExperiment::Compound> filtered_compounds;
    for (const auto& it : compounds)
    {
      // Check if decoy was filtered
      if (std::find(single_decoy_id.begin(), single_decoy_id.end(), it.id) != single_decoy_id.end())
      {
        OPENMS_LOG_DEBUG << "The decoy " << it.id << " was filtered due to missing a respective target." << std::endl;
      }
      else
      {
        filtered_compounds.push_back(it);
      }
    }

    // remove decoy transitions due to missing target
    vector<ReactionMonitoringTransition> filtered_transitions;
    for (const auto& it : exp.getTransitions())
    {
      // Check if compound has any transitions left
      if (std::find(single_decoy_id.begin(), single_decoy_id.end(), it.getCompoundRef()) != single_decoy_id.end())
      {
        OPENMS_LOG_DEBUG << "The decoy " << it.getCompoundRef()
                         << " was filtered due to missing a respective target." << std::endl;
      }
      else
      {
        filtered_transitions.push_back(it);
      }
    }
    exp.setCompounds(filtered_compounds);
    exp.setTransitions(filtered_transitions);
  }

  // =====================================================================
  // Light (memory-efficient) versions of the above methods
  // =====================================================================

  void MRMAssay::reannotateTransitionsLight(OpenSwath::LightTargetedExperiment& exp,
                                            double precursor_mz_threshold,
                                            double product_mz_threshold,
                                            const std::vector<String>& fragment_types,
                                            const std::vector<size_t>& fragment_charges,
                                            bool enable_specific_losses,
                                            bool enable_unspecific_losses,
                                            int round_decPow)
  {
    std::vector<OpenSwath::LightTransition> transitions;
    OpenMS::MRMIonSeries mrmis;

    // Build compound map for quick lookup
    std::map<std::string, const OpenSwath::LightCompound*> compound_map;
    for (const auto& compound : exp.compounds)
    {
      compound_map[compound.id] = &compound;
    }

    // Group transitions by peptide_ref
    std::map<std::string, std::vector<OpenSwath::LightTransition*>> peptide_trans_map;
    for (auto& tr : exp.transitions)
    {
      peptide_trans_map[tr.peptide_ref].push_back(&tr);
    }

    Size progress = 0;
    startProgress(0, exp.transitions.size(), "Annotating transitions (Light)");

    for (auto& pep_it : peptide_trans_map)
    {
      const std::string& peptide_ref = pep_it.first;

      // Get compound
      auto comp_it = compound_map.find(peptide_ref);
      if (comp_it == compound_map.end() || !comp_it->second->isPeptide())
      {
        // Skip non-peptide compounds (metabolomics) - just keep transitions as-is
        for (auto* tr : pep_it.second)
        {
          setProgress(++progress);
          transitions.push_back(*tr);
        }
        continue;
      }

      const OpenSwath::LightCompound* compound = comp_it->second;
      OpenMS::AASequence target_peptide_sequence;
      try
      {
        target_peptide_sequence = AASequence::fromString(compound->sequence);
      }
      catch (Exception::InvalidValue&)
      {
        // Can't parse sequence, keep transitions as-is
        for (auto* tr : pep_it.second)
        {
          setProgress(++progress);
          transitions.push_back(*tr);
        }
        continue;
      }

      int precursor_charge = compound->charge > 0 ? compound->charge : 1;

      MRMIonSeries::IonSeries target_ionseries = mrmis.getIonSeries(
          target_peptide_sequence, precursor_charge, fragment_types,
          fragment_charges, enable_specific_losses,
          enable_unspecific_losses, round_decPow);

      // Generate theoretical precursor m/z
      double precursor_mz = target_peptide_sequence.getMZ(precursor_charge);
      precursor_mz = Math::roundDecimal(precursor_mz, round_decPow);

      for (auto* tr : pep_it.second)
      {
        setProgress(++progress);

        // Annotate transition from theoretical ion series
        std::pair<String, double> targetion = mrmis.annotateIon(target_ionseries, tr->product_mz, product_mz_threshold);

        // Ensure that precursor m/z is within threshold
        if (std::fabs(tr->precursor_mz - precursor_mz) > precursor_mz_threshold)
        {
          targetion.first = "unannotated";
        }

        // Skip unannotated transitions
        if (targetion.first == "unannotated")
        {
          OPENMS_LOG_DEBUG << "[unannotated] Skipping " << target_peptide_sequence.toString()
            << " PrecursorMZ: " << tr->precursor_mz << " ProductMZ: " << tr->product_mz << std::endl;
          continue;
        }

        OpenSwath::LightTransition new_tr = *tr;

        // Set precursor m/z to theoretical value
        new_tr.precursor_mz = precursor_mz;

        // Set product m/z to theoretical value
        new_tr.product_mz = targetion.second;

        // Parse fragment type and number from annotation
        // Format examples: "y7", "b3", "y5-H2O1"
        std::string annotation = targetion.first;
        if (!annotation.empty())
        {
          new_tr.fragment_type = annotation.substr(0, 1);

          // Extract fragment number
          std::string num_str;
          for (size_t i = 1; i < annotation.size(); ++i)
          {
            if (std::isdigit(annotation[i]))
            {
              num_str += annotation[i];
            }
            else
            {
              break;
            }
          }
          if (!num_str.empty())
          {
            new_tr.fragment_nr = std::stoi(num_str);
          }
          new_tr.annotation = annotation;
        }

        transitions.push_back(new_tr);
      }
    }
    endProgress();

    exp.transitions = std::move(transitions);
  }

  void MRMAssay::restrictTransitionsLight(OpenSwath::LightTargetedExperiment& exp,
                                          double lower_mz_limit,
                                          double upper_mz_limit,
                                          const std::vector<std::pair<double, double> >& swathes)
  {
    std::vector<OpenSwath::LightTransition> transitions;

    // Build compound map for quick lookup
    std::map<std::string, const OpenSwath::LightCompound*> compound_map;
    for (const auto& compound : exp.compounds)
    {
      compound_map[compound.id] = &compound;
    }

    Size progress = 0;
    startProgress(0, exp.transitions.size(), "Restricting transitions (Light)");

    for (const auto& tr : exp.transitions)
    {
      setProgress(++progress);

      // Check if product m/z falls into swath from precursor m/z and if yes, skip
      if (!swathes.empty())
      {
        if (MRMAssay::isInSwath_(swathes, tr.precursor_mz, tr.product_mz))
        {
          OPENMS_LOG_DEBUG << "[swath] Skipping PrecursorMZ: " << tr.precursor_mz << " ProductMZ: " << tr.product_mz << std::endl;
          continue;
        }
      }

      // Check if product m/z is outside of m/z boundaries and if yes, skip
      if (tr.product_mz < lower_mz_limit || tr.product_mz > upper_mz_limit)
      {
        OPENMS_LOG_DEBUG << "[mz_limit] Skipping PrecursorMZ: " << tr.precursor_mz << " ProductMZ: " << tr.product_mz << std::endl;
        continue;
      }

      // Append transition
      transitions.push_back(tr);
    }

    exp.transitions = std::move(transitions);
    endProgress();
  }

  void MRMAssay::detectingTransitionsLight(OpenSwath::LightTargetedExperiment& exp,
                                           int min_transitions,
                                           int max_transitions)
  {
    std::vector<OpenSwath::LightTransition> transitions;
    std::vector<OpenSwath::LightCompound> compounds;
    std::vector<OpenSwath::LightProtein> proteins;

    std::unordered_set<std::string> peptide_ids;
    std::unordered_set<std::string> protein_list;

    // Group transitions by peptide_ref
    std::map<std::string, std::vector<OpenSwath::LightTransition>> transitions_map;
    for (const auto& tr : exp.transitions)
    {
      transitions_map[tr.peptide_ref].push_back(tr);
    }

    Size progress = 0;
    startProgress(0, transitions_map.size() + exp.compounds.size() + exp.proteins.size(), "Select detecting transitions (Light)");

    for (auto& m : transitions_map)
    {
      setProgress(++progress);

      // Ensure that all precursors have the minimum number of transitions
      if (m.second.size() >= static_cast<size_t>(min_transitions))
      {
        // LibraryIntensity stores all reference transition intensities of a precursor
        std::vector<double> library_intensity;
        for (const auto& tr : m.second)
        {
          library_intensity.push_back(tr.library_intensity);
        }

        // Sort by intensity, reverse and delete all elements after max_transitions
        std::sort(library_intensity.begin(), library_intensity.end(), std::greater<double>());
        if (static_cast<size_t>(max_transitions) < library_intensity.size())
        {
          library_intensity.resize(max_transitions);
        }

        // Check if transitions are among the ones with maximum intensity
        size_t j = 0;
        for (auto& tr : m.second)
        {
          if (std::find(library_intensity.begin(), library_intensity.end(), tr.library_intensity) != library_intensity.end() &&
              !tr.decoy &&
              j < static_cast<size_t>(max_transitions))
          {
            tr.detecting_transition = true;
            transitions.push_back(tr);
            peptide_ids.insert(tr.peptide_ref);
            j++;
          }
        }
      }
    }

    // Build compound map
    std::map<std::string, const OpenSwath::LightCompound*> compound_map;
    for (const auto& compound : exp.compounds)
    {
      compound_map[compound.id] = &compound;
    }

    for (const auto& compound : exp.compounds)
    {
      setProgress(++progress);

      if (peptide_ids.find(compound.id) != peptide_ids.end())
      {
        compounds.push_back(compound);
        for (const auto& protein_ref : compound.protein_refs)
        {
          protein_list.insert(protein_ref);
        }
      }
    }

    for (const auto& protein : exp.proteins)
    {
      setProgress(++progress);

      if (protein_list.find(protein.id) != protein_list.end())
      {
        proteins.push_back(protein);
      }
    }

    exp.transitions = std::move(transitions);
    exp.compounds = std::move(compounds);
    exp.proteins = std::move(proteins);

    endProgress();
  }

  // =====================================================================
  // Light (memory-efficient) versions of IPF methods
  // =====================================================================

  void MRMAssay::generateTargetInSilicoMapLight_(const OpenSwath::LightTargetedExperiment& exp,
                                                 const std::vector<String>& fragment_types,
                                                 const std::vector<size_t>& fragment_charges,
                                                 bool enable_specific_losses,
                                                 bool enable_unspecific_losses,
                                                 bool enable_ms2_precursors,
                                                 const std::vector<std::pair<double, double> >& swathes,
                                                 int round_decPow,
                                                 size_t max_num_alternative_localizations,
                                                 SequenceMapT& TargetSequenceMap,
                                                 IonMapT& TargetIonMap,
                                                 PeptideMapT& TargetPeptideMap)
  {
    OpenMS::MRMIonSeries mrmis;

    Size progress = 0;
    startProgress(0, exp.compounds.size(), "Generation of target in silico peptide map (Light)");

    for (const auto& compound : exp.compounds)
    {
      setProgress(progress++);

      // Skip non-peptide compounds (metabolites)
      if (!compound.isPeptide())
      {
        continue;
      }

      // Parse sequence - temporary AASequence for peptidoform generation
      OpenMS::AASequence peptide_sequence;
      try
      {
        peptide_sequence = AASequence::fromString(compound.sequence);
      }
      catch (Exception::InvalidValue&)
      {
        OPENMS_LOG_DEBUG << "[uis] Skipping compound (cannot parse sequence): " << compound.id << std::endl;
        continue;
      }

      int precursor_charge = compound.charge > 0 ? compound.charge : 1;
      double precursor_mz = peptide_sequence.getMZ(precursor_charge);
      int precursor_swath = getSwath_(swathes, precursor_mz);

      // Compute all alternative peptidoforms compatible with ModificationsDB
      const std::vector<AASequence> alternative_peptide_sequences = generateTheoreticalPeptidoforms_(peptide_sequence);

      // Some permutations might be too complex, skip if threshold is reached
      if (alternative_peptide_sequences.size() > max_num_alternative_localizations)
      {
        OPENMS_LOG_DEBUG << "[uis] Peptide skipped (too many permutations possible): " << compound.id << std::endl;
        continue;
      }

      // Iterate over all peptidoforms
      for (const auto& alt_aa : alternative_peptide_sequences)
      {
        // Append peptidoform to index
        TargetSequenceMap[precursor_swath][alt_aa.toUnmodifiedString()].insert(alt_aa.toString());

        // Generate theoretical ion series
        auto ionseries = mrmis.getIonSeries(alt_aa, precursor_charge,
            fragment_types, fragment_charges, enable_specific_losses,
            enable_unspecific_losses);

        if (enable_ms2_precursors)
        {
          // Add precursor to theoretical transitions
          double prec_mz = Math::roundDecimal(precursor_mz, round_decPow);
          TargetIonMap[precursor_swath][alt_aa.toUnmodifiedString()].emplace_back(prec_mz, alt_aa.toString());
          TargetPeptideMap[compound.id].emplace_back("MS2_Precursor_i0", prec_mz);
        }

        // Iterate over all theoretical transitions
        for (const auto& im_it : ionseries)
        {
          // Append transition to indices to find interfering transitions
          double fragment_mz = Math::roundDecimal(im_it.second, round_decPow);
          TargetIonMap[precursor_swath][alt_aa.toUnmodifiedString()].emplace_back(fragment_mz, alt_aa.toString());
          TargetPeptideMap[compound.id].emplace_back(im_it.first, fragment_mz);
        }
      }
    }
    endProgress();
  }

  void MRMAssay::generateDecoyInSilicoMapLight_(const OpenSwath::LightTargetedExperiment& exp,
                                                const std::vector<String>& fragment_types,
                                                const std::vector<size_t>& fragment_charges,
                                                bool enable_specific_losses,
                                                bool enable_unspecific_losses,
                                                bool enable_ms2_precursors,
                                                const std::vector<std::pair<double, double> >& swathes,
                                                int round_decPow,
                                                TargetDecoyMapLightT& TargetDecoyMap,
                                                const PeptideMapT& TargetPeptideMap,
                                                const std::map<String, String>& DecoySequenceMap,
                                                IonMapT& DecoyIonMap,
                                                PeptideMapT& DecoyPeptideMap)
  {
    OpenMS::MRMIonSeries mrmis;

    Size progress = 0;
    startProgress(0, exp.compounds.size(), "Generation of decoy in silico peptide map (Light)");

    for (const auto& compound : exp.compounds)
    {
      setProgress(progress++);

      // Skip non-peptide compounds (metabolites)
      if (!compound.isPeptide())
      {
        continue;
      }

      // Skip compounds that are not in target peptide map (e.g., those with too many permutations)
      if (TargetPeptideMap.find(compound.id) == TargetPeptideMap.end())
      {
        continue;
      }

      // Parse sequence
      OpenMS::AASequence peptide_sequence;
      try
      {
        peptide_sequence = AASequence::fromString(compound.sequence);
      }
      catch (Exception::InvalidValue&)
      {
        continue;
      }

      String unmodified_sequence = peptide_sequence.toUnmodifiedString();

      // Check if we have a decoy sequence for this target
      auto decoy_it = DecoySequenceMap.find(unmodified_sequence);
      if (decoy_it == DecoySequenceMap.end())
      {
        continue;
      }

      // Create decoy compound
      OpenSwath::LightCompound decoy_compound = compound;
      decoy_compound.id = "DECOY_" + compound.id;

      // Build decoy sequence with modifications
      AASequence decoy_sequence = AASequence::fromString(decoy_it->second);

      // Generate theoretical peptidoforms for decoy using target as template
      const std::vector<AASequence> alternative_peptide_sequences =
          generateTheoreticalPeptidoformsDecoy_(peptide_sequence, decoy_sequence);

      int precursor_charge = compound.charge > 0 ? compound.charge : 1;

      // Use TARGET peptide's precursor for SWATH index (same as heavy version)
      double target_precursor_mz = peptide_sequence.getMZ(precursor_charge);
      int precursor_swath = getSwath_(swathes, target_precursor_mz);

      // Iterate over all peptidoforms
      for (const auto& alt_aa : alternative_peptide_sequences)
      {
        // Generate theoretical ion series
        auto ionseries = mrmis.getIonSeries(alt_aa, precursor_charge,
            fragment_types, fragment_charges, enable_specific_losses,
            enable_unspecific_losses);

        if (enable_ms2_precursors)
        {
          // Use TARGET precursor m/z (same as heavy version)
          double prec_mz = Math::roundDecimal(target_precursor_mz, round_decPow);
          DecoyIonMap[precursor_swath][alt_aa.toUnmodifiedString()].emplace_back(prec_mz, alt_aa.toString());
          // Use TARGET compound id as key (same as heavy version)
          DecoyPeptideMap[compound.id].emplace_back("MS2_Precursor_i0", prec_mz);
        }

        for (const auto& im_it : ionseries)
        {
          double fragment_mz = Math::roundDecimal(im_it.second, round_decPow);
          DecoyIonMap[precursor_swath][alt_aa.toUnmodifiedString()].emplace_back(fragment_mz, alt_aa.toString());
          // Use TARGET compound id as key (same as heavy version)
          DecoyPeptideMap[compound.id].emplace_back(im_it.first, fragment_mz);
        }
      }

      // Store decoy compound in map
      decoy_compound.sequence = alternative_peptide_sequences.empty() ? decoy_it->second : alternative_peptide_sequences[0].toString();
      TargetDecoyMap[compound.id] = decoy_compound;
    }
    endProgress();
  }

  void MRMAssay::generateTargetAssaysLight_(const OpenSwath::LightTargetedExperiment& exp,
                                            std::vector<OpenSwath::LightTransition>& transitions,
                                            double mz_threshold,
                                            const std::vector<std::pair<double, double> >& swathes,
                                            int round_decPow,
                                            const PeptideMapT& TargetPeptideMap,
                                            const IonMapT& TargetIonMap)
  {
    // Build compound lookup map
    std::map<std::string, const OpenSwath::LightCompound*> compound_map;
    for (const auto& compound : exp.compounds)
    {
      compound_map[compound.id] = &compound;
    }

    Size progress = 0;
    startProgress(0, TargetPeptideMap.size(), "Generation of target identification transitions (Light)");

    size_t transition_idx = transitions.size();

    for (const auto& pep_it : TargetPeptideMap)
    {
      setProgress(progress++);

      const String& peptide_id = pep_it.first;
      auto comp_it = compound_map.find(peptide_id);
      if (comp_it == compound_map.end())
      {
        continue;
      }

      const OpenSwath::LightCompound* compound = comp_it->second;

      OpenMS::AASequence peptide_sequence;
      try
      {
        peptide_sequence = AASequence::fromString(compound->sequence);
      }
      catch (Exception::InvalidValue&)
      {
        continue;
      }

      int precursor_charge = compound->charge > 0 ? compound->charge : 1;
      double precursor_mz = peptide_sequence.getMZ(precursor_charge);
      int precursor_swath = getSwath_(swathes, precursor_mz);
      String unmod_str = peptide_sequence.toUnmodifiedString();

      // Check if we have data for this SWATH/sequence combination
      auto swath_it = TargetIonMap.find(precursor_swath);
      if (swath_it == TargetIonMap.end())
      {
        continue;
      }
      auto seq_it = swath_it->second.find(unmod_str);
      if (seq_it == swath_it->second.end())
      {
        continue;
      }

      // Sort and deduplicate transitions
      IonSeries sorted_ions = pep_it.second;
      std::sort(sorted_ions.begin(), sorted_ions.end());
      sorted_ions.erase(std::unique(sorted_ions.begin(), sorted_ions.end()), sorted_ions.end());

      for (const auto& ion : sorted_ions)
      {
        const std::string& annotation = ion.first;
        double product_mz = ion.second;

        // Find matching peptidoforms
        std::vector<std::string> matching_pfs = getMatchingPeptidoforms_(
            product_mz, seq_it->second, mz_threshold);

        if (matching_pfs.empty())
        {
          continue;
        }

        // Create light transition
        OpenSwath::LightTransition tr;
        tr.transition_name = String(transition_idx++) + "_UIS_" + peptide_id + "_" +
            String(Math::roundDecimal(precursor_mz, round_decPow)) + "_" +
            String(Math::roundDecimal(product_mz, round_decPow));
        tr.peptide_ref = peptide_id;
        tr.precursor_mz = precursor_mz;
        tr.product_mz = product_mz;
        tr.library_intensity = 1.0;
        tr.decoy = false;
        tr.detecting_transition = false;
        tr.identifying_transition = true;
        tr.quantifying_transition = false;
        tr.annotation = annotation;
        tr.peptidoforms = matching_pfs;

        // Parse fragment type/number from annotation
        if (!annotation.empty() && annotation != "MS2_Precursor_i0")
        {
          tr.fragment_type = annotation.substr(0, 1);
          std::string num_str;
          for (size_t i = 1; i < annotation.size() && std::isdigit(annotation[i]); ++i)
          {
            num_str += annotation[i];
          }
          if (!num_str.empty())
          {
            tr.fragment_nr = std::stoi(num_str);
          }
        }

        transitions.push_back(std::move(tr));
      }
    }
    endProgress();
  }

  void MRMAssay::generateDecoyAssaysLight_(const OpenSwath::LightTargetedExperiment& exp,
                                           std::vector<OpenSwath::LightTransition>& transitions,
                                           double mz_threshold,
                                           const std::vector<std::pair<double, double> >& swathes,
                                           int round_decPow,
                                           const PeptideMapT& DecoyPeptideMap,
                                           const TargetDecoyMapLightT& TargetDecoyMap,
                                           const IonMapT& DecoyIonMap,
                                           const IonMapT& TargetIonMap)
  {
    // Build compound lookup map for targets
    std::map<std::string, const OpenSwath::LightCompound*> compound_map;
    for (const auto& compound : exp.compounds)
    {
      compound_map[compound.id] = &compound;
    }

    Size progress = 0;
    startProgress(0, DecoyPeptideMap.size(), "Generation of decoy identification transitions (Light)");

    size_t transition_idx = transitions.size();

    for (const auto& pep_it : DecoyPeptideMap)
    {
      setProgress(progress++);

      // DecoyPeptideMap keys are TARGET compound ids (same as heavy version)
      const String& target_id = pep_it.first;

      auto decoy_comp_it = TargetDecoyMap.find(target_id);
      if (decoy_comp_it == TargetDecoyMap.end())
      {
        continue;
      }

      const OpenSwath::LightCompound& decoy_compound = decoy_comp_it->second;

      auto target_comp_it = compound_map.find(target_id);
      if (target_comp_it == compound_map.end())
      {
        continue;
      }
      const OpenSwath::LightCompound* target_compound = target_comp_it->second;

      OpenMS::AASequence decoy_sequence;
      OpenMS::AASequence target_sequence;
      try
      {
        decoy_sequence = AASequence::fromString(decoy_compound.sequence);
        target_sequence = AASequence::fromString(target_compound->sequence);
      }
      catch (Exception::InvalidValue&)
      {
        continue;
      }

      int precursor_charge = decoy_compound.charge > 0 ? decoy_compound.charge : 1;
      double decoy_precursor_mz = decoy_sequence.getMZ(precursor_charge);
      double target_precursor_mz = target_sequence.getMZ(precursor_charge);

      int decoy_precursor_swath = getSwath_(swathes, decoy_precursor_mz);
      int target_precursor_swath = getSwath_(swathes, target_precursor_mz);

      String decoy_unmod_str = decoy_sequence.toUnmodifiedString();
      String target_unmod_str = target_sequence.toUnmodifiedString();

      // Check if we have decoy ion data
      auto decoy_swath_it = DecoyIonMap.find(decoy_precursor_swath);
      if (decoy_swath_it == DecoyIonMap.end())
      {
        continue;
      }
      auto decoy_seq_it = decoy_swath_it->second.find(decoy_unmod_str);
      if (decoy_seq_it == decoy_swath_it->second.end())
      {
        continue;
      }

      // Sort and deduplicate transitions
      IonSeries sorted_ions = pep_it.second;
      std::sort(sorted_ions.begin(), sorted_ions.end());
      sorted_ions.erase(std::unique(sorted_ions.begin(), sorted_ions.end()), sorted_ions.end());

      for (const auto& ion : sorted_ions)
      {
        const std::string& annotation = ion.first;
        double product_mz = ion.second;

        // Find matching decoy peptidoforms
        std::vector<std::string> matching_pfs = getMatchingPeptidoforms_(
            product_mz, decoy_seq_it->second, mz_threshold);

        if (matching_pfs.empty())
        {
          continue;
        }

        // Check for overlap with target transitions
        auto target_swath_it = TargetIonMap.find(target_precursor_swath);
        if (target_swath_it != TargetIonMap.end())
        {
          auto target_seq_it = target_swath_it->second.find(target_unmod_str);
          if (target_seq_it != target_swath_it->second.end())
          {
            std::vector<std::string> target_overlap = getMatchingPeptidoforms_(
                product_mz, target_seq_it->second, mz_threshold);
            if (!target_overlap.empty())
            {
              // Skip this decoy transition as it overlaps with target
              continue;
            }
          }
        }

        // Create light transition
        OpenSwath::LightTransition tr;
        tr.transition_name = String(transition_idx++) + "_UISDECOY_" + decoy_compound.id + "_" +
            String(Math::roundDecimal(decoy_precursor_mz, round_decPow)) + "_" +
            String(Math::roundDecimal(product_mz, round_decPow));
        tr.peptide_ref = decoy_compound.id;
        tr.precursor_mz = decoy_precursor_mz;
        tr.product_mz = product_mz;
        tr.library_intensity = 1.0;
        tr.decoy = true;
        tr.detecting_transition = false;
        tr.identifying_transition = true;
        tr.quantifying_transition = false;
        tr.annotation = annotation;
        tr.peptidoforms = matching_pfs;

        // Parse fragment type/number from annotation
        if (!annotation.empty() && annotation != "MS2_Precursor_i0")
        {
          tr.fragment_type = annotation.substr(0, 1);
          std::string num_str;
          for (size_t i = 1; i < annotation.size() && std::isdigit(annotation[i]); ++i)
          {
            num_str += annotation[i];
          }
          if (!num_str.empty())
          {
            tr.fragment_nr = std::stoi(num_str);
          }
        }

        transitions.push_back(std::move(tr));
      }
    }
    endProgress();
  }

  void MRMAssay::uisTransitionsLight(OpenSwath::LightTargetedExperiment& exp,
                                     const std::vector<String>& fragment_types,
                                     const std::vector<size_t>& fragment_charges,
                                     bool enable_specific_losses,
                                     bool enable_unspecific_losses,
                                     bool enable_ms2_precursors,
                                     double mz_threshold,
                                     const std::vector<std::pair<double, double> >& swathes,
                                     int round_decPow,
                                     size_t max_num_alternative_localizations,
                                     int shuffle_seed,
                                     bool disable_decoy_transitions)
  {
    // Use same map types as heavy version - they use std::string keys
    SequenceMapT TargetSequenceMap;
    IonMapT TargetIonMap;
    PeptideMapT TargetPeptideMap;

    std::vector<OpenSwath::LightTransition> transitions;

    // Step 1: Generate target in-silico map
    generateTargetInSilicoMapLight_(exp, fragment_types, fragment_charges,
        enable_specific_losses, enable_unspecific_losses, enable_ms2_precursors,
        swathes, round_decPow, max_num_alternative_localizations,
        TargetSequenceMap, TargetIonMap, TargetPeptideMap);

    // Step 2: Generate target identification transitions
    generateTargetAssaysLight_(exp, transitions, mz_threshold, swathes,
        round_decPow, TargetPeptideMap, TargetIonMap);

    if (!disable_decoy_transitions)
    {
      // Step 3: Generate decoy sequences and in-silico map
      std::map<String, String> DecoySequenceMap;
      IonMapT DecoyIonMap;
      PeptideMapT DecoyPeptideMap;
      TargetDecoyMapLightT TargetDecoyMap;

      // Reuse the existing heavy generateDecoySequences_ - it only uses SequenceMapT
      generateDecoySequences_(TargetSequenceMap, DecoySequenceMap, shuffle_seed);

      generateDecoyInSilicoMapLight_(exp, fragment_types, fragment_charges,
          enable_specific_losses, enable_unspecific_losses, enable_ms2_precursors,
          swathes, round_decPow, TargetDecoyMap, TargetPeptideMap,
          DecoySequenceMap, DecoyIonMap, DecoyPeptideMap);

      // Step 4: Generate decoy identification transitions
      generateDecoyAssaysLight_(exp, transitions, mz_threshold, swathes,
          round_decPow, DecoyPeptideMap, TargetDecoyMap, DecoyIonMap, TargetIonMap);
    }

    // Append new transitions to experiment
    for (auto& tr : transitions)
    {
      exp.transitions.push_back(std::move(tr));
    }
  }

}

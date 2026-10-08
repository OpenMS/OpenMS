// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Lucia Espona, Moritz Freidank $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------



namespace OpenMS
{

/**
    @brief Resolves ambiguous annotations of features with peptide identifications.

    The peptide identifications are filtered so that only one identification
    with a single hit (with the best score) is associated to each feature.
    (If two IDs have the same best score, the first one is selected.)

    The map functions work on the identification data of a map and the features' links to it: an identification
    that a feature no longer links is unassigned, and a match removed from an identification is removed from the
    identification data. Scores are the primary scores of the runs. Peptide identifications of a map are converted
    to identification data first and back afterwards (IdentificationDataConverter), so they are resolved the same
    way; their order follows the identification data (the order of the identifications when they were converted).
*/
class OPENMS_DLLAPI IDConflictResolverAlgorithm
{
public:
  /** @brief Resolves ambiguous annotations of features with peptide identifications.
    The the filtered identifications are added to the vector of unassigned peptides
    and also reduced to a single best hit.

    @param[in] features Features to work on
    @param[in,out] keep_matching Keeps all IDs that match the modified sequence of the best
    hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
  **/
  static void resolve(FeatureMap& features, bool keep_matching = false);

  /** @brief Resolves ambiguous annotations of consensus features with peptide identifications.
    The the filtered identifications are added to the vector of unassigned peptides
    and also reduced to a single best hit.
    
    @param[in] features Features to work on
    @param[in,out] keep_matching Keeps all IDs that match the modified sequence of the best
    hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
  **/
  static void resolve(ConsensusMap& features, bool keep_matching = false);

  /** @brief Resolves ambiguous annotations of features with peptide identifications using rank aggregation.

    For each feature, peptide hits across all identifications are aggregated by rank.
    Each unique sequence is assigned a rank in every identification in which it appears
    (rank 0 = best hit, 1 = second best, etc.). Sequences not found in an identification
    receive a penalty rank equal to the maximum number of considered hits. The aggregate
    score for each sequence is computed as:
    @code
      1.0 - (sum_of_ranks + penalty_for_missing_runs) / (max_hits * n_runs)
    @endcode
    The sequence with the highest aggregate score is selected as the winner and the
    corresponding best-scoring identification is kept with only that hit.
    All other identifications are moved to the unassigned list.

    @param[in] features FeatureMap to work on
  **/
  static void resolveAllHitRankAggregation(FeatureMap& features);

  /** @brief Resolves ambiguous annotations of consensus features with peptide identifications using rank aggregation.

    For each consensus feature, peptide hits across all identifications are aggregated by rank.
    Each unique sequence is assigned a rank in every identification in which it appears
    (rank 0 = best hit, 1 = second best, etc.). Sequences not found in an identification
    receive a penalty rank equal to the maximum number of considered hits. The aggregate
    score for each sequence is computed as:
    @code
      1.0 - (sum_of_ranks + penalty_for_missing_runs) / (max_hits * n_runs)
    @endcode
    The sequence with the highest aggregate score is selected as the winner and the
    corresponding best-scoring identification is kept with only that hit.
    All other identifications are moved to the unassigned list.

    @param[in] features ConsensusMap to work on
  **/
  static void resolveAllHitRankAggregation(ConsensusMap& features);

  /** @brief In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear.
   This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
   combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
   appears only once, i.e. no multiplicities.
   **/
  static void resolveBetweenFeatures(FeatureMap& features);
  
  /** @brief In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear.
   This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
   combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
   appears only once, i.e. no multiplicities.
   **/
  static void resolveBetweenFeatures(ConsensusMap& features);

  /// What reduceToOnePerSpectrum() found and did. All counts are over the input list.
  struct OPENMS_DLLAPI UnresolvedIdentifications
  {
    /// Identifications removed because another one claimed the same
    /// (spectrum reference, top-hit sequence, top-hit charge).
    Size removed = 0;
    /// Spectra that still carry more than one identification AFTER the reduction, i.e. ones
    /// whose identifications name different peptidoforms (a chimeric spectrum, or search
    /// engines that disagree). Not reduced - both are quantifiable - but worth reporting.
    Size multiply_identified_spectra = 0;
    /// Identifications that could not be keyed because they carry no spectrum reference.
    /// Left untouched: without a reference there is nothing to tell them apart by.
    Size without_spectrum_reference = 0;
    /// Identifications whose group disagreed on isHigherScoreBetter(), so "best" was undefined.
    /// Left untouched rather than reduced by a coin flip.
    Size inconsistent_score_direction = 0;
    /// First reduced group, as "<spectrum reference> / <sequence> / charge <n>", for logging.
    std::string example;
  };

  /**
    @brief Reduce identifications a quantification workflow cannot tell apart to one per spectrum.

    Quantification tools measure one value per (spectrum, peptidoform, charge). Where the input
    carries several identifications of that same triple - e.g. two search engines that agree,
    concatenated rather than combined - only one of them can own the measurement, and keeping
    all of them counts one measurement several times.

    This keeps the best-scoring identification of each such group and removes the rest. That is a
    tie-break on scores something upstream already assigned, not a scoring decision: combining
    scores from different engines is ConsensusIDAlgorithm's job and is deliberately not done here.

    Identifications are keyed on the top hit as stored; hits are never re-sorted, so the key
    matches what the exporters and quantifiers read. Surviving identifications keep their
    original relative order.

    @param[in,out] ids One identification run's peptide identifications.
    @return What was found and removed; see UnresolvedIdentifications.
  **/
  static UnresolvedIdentifications reduceToOnePerSpectrum(PeptideIdentificationList& ids);

};

}// namespace OpenMS


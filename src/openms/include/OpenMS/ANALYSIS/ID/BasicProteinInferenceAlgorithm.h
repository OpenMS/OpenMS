// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>

namespace OpenMS
{

  /** \brief Algorithm class that implements simple protein inference by aggregation of peptide scores.
   * It has multiple parameter options like the aggregation method, when to distinguish peptidoforms,
   * and if you want to use shared peptides ("use_shared_peptides").
   * First, the best PSM per spectrum is used, then only the best PSM per peptidoform is aggregated.
   * Peptidoforms can optionally be distinguished via the treat_X_separate parameters:
   * - Modifications (modified sequence string)
   * - Charge states
   * The algorithm assumes posteriors or posterior error probabilities and converts to posteriors initially.
   * Possible aggregation methods that can be set via the parameter "aggregation_method" are:
   * - "maximum" (default)
   * - "sum"
   * - "product" (ignoring zeroes)
   * Annotation of the number of peptides used for aggregation can be disabled (see parameters).
   * Supports multiple runs but goes through them one by one iterating over the full PeptideIdentification vector.
   */
  class OPENMS_DLLAPI BasicProteinInferenceAlgorithm :
    public DefaultParamHandler,
    public ProgressLogger
  {
    public:

    typedef std::unordered_map<std::string, std::map<Int, PeptideHit*>> SequenceToChargeToPSM;

    /**
     * @brief The aggregation method
     */
    enum class AggregationMethod
    {
      PROD, ///< aggregate by product (ignore zeroes)
      SUM, ///< aggregate by summing
      BEST ///< aggregate by maximum/minimum
    };

    /// Default constructor
    BasicProteinInferenceAlgorithm();

    /**
     * Performs the actual inference based on best psm per peptide in @p pep_ids per run in @p prot_ids.
     * Sorts and filters psms in @p pep_ids. Annotates results in @p prot_ids.
     * Associations (via getIdentifier) for peptides to protein runs need to be correct.
     */
    void run(std::vector<PeptideIdentification>& pep_ids, std::vector<ProteinIdentification>& prot_ids) const;

    /**
     * Performs the actual inference based on best psm per peptide in @p pep_ids per run in @p prot_id.
     * Sorts and filters psms in @p pep_ids. Annotates results in @p prot_id.
     * Associations (via getIdentifier) for peptides to protein runs need to be correct.
     */
    void run(std::vector<PeptideIdentification>& pep_ids, ProteinIdentification& prot_id) const;

    /**
     * Performs the actual inference based on best psm per peptide in @p cmap for proteins from @p prot_id.
     * Ideally @p prot_id is the union of proteins in all runs of @p cmap.
     * Sorts and filters psms in @p pep_ids. Annotates results in @p prot_id.
     * Associations (via getIdentifier) for peptides to protein runs ARE IGNORED and all pep_ids used.
     * @todo allow checking matching IDs
     */
    void run(ConsensusMap& cmap, ProteinIdentification& prot_id, bool include_unassigned) const;

  private:

    /**
     * @brief Performs simple aggregation-based inference on one protein run.
     * @param acc_to_protein_hitP_and_count Maps Accessions to a pair of ProteinHit pointers
     *  and number of peptidoforms encountered
     * @param best_pep Maps (un)modified peptide sequence to a map from charge (0 when unconsidered) to the
     *  best PeptideHit pointer
     * @param prot_run The current run to process
     * @param pep_ids Peptides for the current run to process
     */
    void processRun_(
      std::unordered_map<std::string, std::pair<ProteinHit*, Size>>& acc_to_protein_hitP_and_count,
      SequenceToChargeToPSM& best_pep,
      ProteinIdentification& prot_run,
      std::vector<PeptideIdentification>& pep_ids) const;

    /**
     * @brief fills and updates the map of best peptide scores @p best_pep (by sequence or modified sequence, depending on algorithm settings)
     * @param best_pep (mod.) sequence to charge to pointer of best PSM (PeptideHit*)
     * @param pep_ids the spectra with PSMs
     * @param overall_score_type the pre-determined type name to raise an error if mixed types occur
     * @param higher_better if for this score type higher is better
     * @param run_id only process peptides associated with this run_id (e.g. proteinID run getIdentifier())
     */
    void aggregatePeptideScores_(
        SequenceToChargeToPSM& best_pep,
        std::vector<PeptideIdentification>& pep_ids,
        const String& overall_score_type,
        bool higher_better,
        const std::string& run_id) const;

    /**
     * @brief aggregates and updates protein scores based on aggregation settings and aggregated peptide level results in
     * prefilled @p best_pep
     * @param acc_to_protein_hitP_and_count the results to fill
     * @param best_pep best psm per peptide to read the score
     * @param pep_scores if the score is a posterior error probability -> Auto-converts to posterior probability
     * @param higher_better if for the score higher is better. Assume score is unconverted.
     */
    void updateProteinScores_(
        std::unordered_map<std::string, std::pair<ProteinHit*, Size>>& acc_to_protein_hitP_and_count,
        const SequenceToChargeToPSM& best_pep,
        bool pep_scores,
        bool higher_better) const;

    /// get the AggregationMethod enum from a @p method_string
    AggregationMethod aggFromString_(const std::string& method_string) const;

    /// check if a @p score_name is compatible to the chosen @p aggregation_method
    /// I.e. only probabilities can be used for multiplication
    void checkCompat_(
        const String& score_type,
        const AggregationMethod& aggregation_method
        ) const;

    /// check if a @p score_type is compatible to the chosen @p aggregation_method
    /// I.e. only probabilities can be used for multiplication
    void checkCompat_(
        const IDScoreSwitcherAlgorithm::ScoreType& score_type,
        const AggregationMethod& aggregation_method
        ) const;

    /// get the initial score value based on the chosen @p aggregation_method, @p higher_better is needed for "best" score
    double getInitScoreForAggMethod_(const AggregationMethod& aggregation_method, bool higher_better) const;

    /// get lambda function to aggregate scores
    typedef double (*fptr)(double, double);
    fptr aggFunFromEnum_(const BasicProteinInferenceAlgorithm::AggregationMethod& agg_method, bool higher_better) const;
  };
} //namespace OpenMS

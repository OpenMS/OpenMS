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
#include <OpenMS/METADATA/PeptideIdentificationList.h>

namespace OpenMS
{

  /**
    @brief Simple protein inference by aggregation of per-peptide PSM scores.

    First takes the best PSM per spectrum, then keeps the best PSM per peptidoform
    (where "peptidoform" is widened or narrowed by the @c treat_charge_variants_separately
    and @c treat_modification_variants_separately parameters), and finally aggregates the
    peptide-level scores onto the proteins using one of the methods exposed via the
    @c score_aggregation_method parameter.

    Configurable behaviour is exposed through @ref DefaultParamHandler — see the defaults
    installed by the constructor for the full list of supported keys, in particular:

      - @c "score_aggregation_method" — one of @c "best", @c "product", @c "sum",
        @c "maximum"; maps onto @ref AggregationMethod via @ref aggFromString_. The
        @c "best" / @c "maximum" string both produce the @c BEST mode.
      - @c "treat_charge_variants_separately" — distinguish charge variants of the same
        modified sequence as distinct peptidoforms (default @c "true").
      - @c "treat_modification_variants_separately" — distinguish modified vs. unmodified
        variants of the same backbone sequence (default @c "true").
      - @c "use_shared_peptides" — if @c "true", shared peptides count as evidence for
        every protein they map to (default @c "true").
      - @c "skip_count_annotation" — if set, the per-peptide count annotation on the
        protein hits is suppressed (default @c "false").
      - @c "annotate_indistinguishable_groups" — compute and annotate indistinguishable
        protein groups (default @c "true").
      - @c "greedy_group_resolution" — resolve shared peptides to a single best protein
        (razor-peptide style; default @c "false").
      - @c "min_peptides_per_protein" — minimum peptide count required for a protein to
        be reported (default 1).
      - @c "score_type" — explicit PSM score type to use; empty falls back to the main
        score.

    The algorithm assumes posteriors or posterior error probabilities; PEPs
    are converted to posteriors as part of scoring. Multiple runs are
    supported, each processed independently.

    @ingroup Analysis_ID
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
      @brief Run inference per protein-ID run, iterating each @p prot_ids entry separately.

      For every entry in @p prot_ids, only peptides whose @c getIdentifier matches that
      run's @c getIdentifier are processed (other peptides are ignored for that run).
      @p pep_ids is sorted and filtered to best-PSM-per-peptidoform; @p prot_ids is
      annotated with the aggregated scores and (unless @c "skip_count_annotation" is set)
      with per-protein peptide counts.

      @param[in,out] pep_ids   Peptide identifications across all runs; sorted/filtered in place.
      @param[in,out] prot_ids  One protein-identification run per entry; scores and per-peptide counts annotated in place.
      @throws Exception::InvalidParameter If PSMs of a peptide carry different score types (mixed score types are not supported).
    */
    void run(PeptideIdentificationList& pep_ids, std::vector<ProteinIdentification>& prot_ids) const;

    /**
      @brief Run inference for a single protein-ID run.

      Convenience overload of the multi-run version: only peptides whose @c getIdentifier
      matches @p prot_id.getIdentifier() are processed; the others are ignored.

      @param[in,out] pep_ids   Peptide identifications for this run; sorted/filtered in place.
      @param[in,out] prot_id   Protein-identification run to annotate with aggregated scores.
      @throws Exception::InvalidParameter If PSMs of a peptide carry different score types.
    */
    void run(PeptideIdentificationList& pep_ids, ProteinIdentification& prot_id) const;

    /**
      @brief Run inference over a @ref ConsensusMap, treating every peptide identification it carries as evidence for the proteins in @p prot_id.

      Differs from the per-run overloads above by **ignoring** the @c getIdentifier
      association between peptides and protein runs — every peptide id in @p cmap (and
      optionally in the unassigned list, see @p include_unassigned) is used. @p prot_id is
      expected to be the union of the proteins of all runs in @p cmap.

      @param[in,out] cmap                Consensus map providing the peptide identifications; PSMs may be sorted/filtered in place.
      @param[in,out] prot_id             Protein-identification run to annotate with aggregated scores.
      @param[in]     include_unassigned  If true, also include @ref ConsensusMap::getUnassignedPeptideIdentifications.
      @throws Exception::InvalidParameter If PSMs of a peptide carry different score types.
      @todo JuliaP Allow checking that peptide / protein IDs reference the same run identifier.
    */
    void run(ConsensusMap& cmap, ProteinIdentification& prot_id, bool include_unassigned) const;

  private:

    /**
     * @brief Performs simple aggregation-based inference on one protein run.
     * @param[in,out] acc_to_protein_hitP_and_count Maps Accessions to a pair of ProteinHit pointers
     *  and number of peptidoforms encountered
     * @param[in,out] best_pep Maps (un)modified peptide sequence to a map from charge (0 when unconsidered) to the
     *  best PeptideHit pointer
     * @param[in,out] prot_run The current run to process
     * @param[in,out] pep_ids Peptides for the current run to process
     */
    void processRun_(
      std::unordered_map<std::string, std::pair<ProteinHit*, Size>>& acc_to_protein_hitP_and_count,
      SequenceToChargeToPSM& best_pep,
      ProteinIdentification& prot_run,
      PeptideIdentificationList& pep_ids) const;

    /**
     * @brief fills and updates the map of best peptide scores @p best_pep (by sequence or modified sequence, depending on algorithm settings)
     * @param[in,out] best_pep (mod.) sequence to charge to pointer of best PSM (PeptideHit*)
     * @param[in,out] pep_ids the spectra with PSMs
     * @param[in] overall_score_type the pre-determined type name to raise an error if mixed types occur
     * @param[in] higher_better if for this score type higher is better
     * @param[in] run_id only process peptides associated with this run_id (e.g. proteinID run getIdentifier())
     */
    void aggregatePeptideScores_(
        SequenceToChargeToPSM& best_pep,
        PeptideIdentificationList& pep_ids,
        const String& overall_score_type,
        bool higher_better,
        const std::string& run_id) const;

    /**
     * @brief aggregates and updates protein scores based on aggregation settings and aggregated peptide level results in
     * prefilled @p best_pep
     * @param[in,out] acc_to_protein_hitP_and_count the results to fill
     * @param[in] best_pep best psm per peptide to read the score
     * @param[in] pep_scores if the score is a posterior error probability -> Auto-converts to posterior probability
     * @param[in] higher_better if for the score higher is better. Assume score is unconverted.
     */
    void updateProteinScores_(
        std::unordered_map<std::string, std::pair<ProteinHit*, Size>>& acc_to_protein_hitP_and_count,
        const SequenceToChargeToPSM& best_pep,
        bool pep_scores,
        bool higher_better) const;

    /**
      @brief Map a @c score_aggregation_method parameter string to the @ref AggregationMethod enum.

      Recognised values: @c "product", @c "sum", @c "best", @c "maximum"
      (@c "best" and @c "maximum" both produce @ref AggregationMethod::BEST).

      @param[in] method_string Parameter string.
      @return Matching @ref AggregationMethod.
    */
    AggregationMethod aggFromString_(const std::string& method_string) const;

    /**
      @brief Reject score-type / aggregation-method combinations that don't make statistical sense.

      Multiplication (@ref AggregationMethod::PROD) is only meaningful for probability-typed
      scores; other combinations either throw or log a warning. Uses the score-type @em name.

      @param[in] score_type        Name of the PSM score type.
      @param[in] aggregation_method Aggregation mode chosen via @c "score_aggregation_method".
    */
    void checkCompat_(
        const String& score_type,
        const AggregationMethod& aggregation_method
        ) const;

    /**
      @brief Same as the string overload, but takes a typed @ref IDScoreSwitcherAlgorithm::ScoreType so the check can be done after the score-switcher has classified the score.

      @param[in] score_type        Typed score classification.
      @param[in] aggregation_method Aggregation mode chosen via @c "score_aggregation_method".
    */
    void checkCompat_(
        const IDScoreSwitcherAlgorithm::ScoreType& score_type,
        const AggregationMethod& aggregation_method
        ) const;

    /**
      @brief Return the identity-element initial score for the chosen aggregation method.

      For example, 0 for SUM, 1 for PROD, and the worst-possible score (depending on
      @p higher_better) for BEST.

      @param[in] aggregation_method Aggregation mode chosen via @c "score_aggregation_method".
      @param[in] higher_better      Whether higher score values are better (used only for BEST).
      @return Initial accumulator value to start the aggregation from.
    */
    double getInitScoreForAggMethod_(const AggregationMethod& aggregation_method, bool higher_better) const;

    /// Function-pointer type for a two-argument score accumulator
    typedef double (*fptr)(double, double);

    /**
      @brief Pick the two-argument accumulator function matching the chosen aggregation method.

      @param[in] agg_method   Aggregation mode chosen via @c "score_aggregation_method".
      @param[in] higher_better Whether higher score values are better (used only for BEST, to pick max vs min).
      @return Function pointer with signature @c double(double,double).
    */
    fptr aggFunFromEnum_(const BasicProteinInferenceAlgorithm::AggregationMethod& agg_method, bool higher_better) const;
  };
} //namespace OpenMS

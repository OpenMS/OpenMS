// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <unordered_map>

#include <vector>
#include <unordered_set>

namespace OpenMS
{

  struct ScoreToTgtDecLabelPairs;

  /**
    @brief Calculates false discovery rates (FDR) from identifications

    Either two runs of forward and decoy database identification or one run containing both (with annotations) can be used to annotate each of the peptide hits with an FDR or q-value.

    q-values are basically only adjusted p-values, also ranging from 0 to 1, with lower values being preferable.
    When looking at the list of hits ordered by q-values, then a specific q-value of @em x means that @em x*100 percent of hits with a q-value <= @em x are expected to be false positives.

    Only simple target-decoy FDRs are supported with a formula depending on the "conservative" parameter:
    - false: (D+1)/T.
    - true: (D+1)/(T+D) [for comparison with protein level FDR used by other tools like e.g., Fido]
    For protein groups, a group is considered as a target when it contains at least one target protein.
    Group level FDRs assume the same score type as on protein level.

    For peptide hits, a hit is considered target also if it maps to both
    a target and a decoy protein (i.e. "target+decoy") as value in the
    "target_decoy" metavalue e.g. annotated by @ref TOPP_PeptideIndexer

    @note The parameter add_decoy_proteins currently does not affect groups

    @htmlinclude OpenMS_FalseDiscoveryRate.parameters

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI FalseDiscoveryRate :
    public DefaultParamHandler
  {
public:
    ///Default constructor
    FalseDiscoveryRate();

    /**
       @brief Calculates the FDR of two runs, a forward run and a decoy run on peptide level

       @param fwd_ids forward peptide identifications
       @param rev_ids reverse peptide identifications
    */
    void apply(std::vector<PeptideIdentification>& fwd_ids, std::vector<PeptideIdentification>& rev_ids) const;

    /**
    @brief Calculates the FDR of one run from a concatenated sequence DB search.    

    @param id peptide identifications, containing target and decoy hits
    @param annotate_peptide_fdr adds the peptide q-value or peptide fdr meta value to each PSM. Calculation uses best PSM per peptide.
    */
    void apply(std::vector<PeptideIdentification>& id, bool annotate_peptide_fdr = false) const;

    /**
    @brief Calculates the FDR of two runs, a forward run and decoy run on protein level

    @param fwd_ids forward protein identifications
    @param rev_ids reverse protein identifications
    */
    void apply(std::vector<ProteinIdentification>& fwd_ids, std::vector<ProteinIdentification>& rev_ids) const;

    /**
    @brief Calculate the FDR of one run from a concatenated sequence db search

    @param ids protein identifications, containing target and decoy hits
    */
    void apply(std::vector<ProteinIdentification>& ids) const;

    /**
    @brief Calculate the FDR based on PEPs or PPs (if present) and modifies the IDs inplace

    @param ids protein identifications, containing PEP scores (not necessarily) annotated with target decoy.
    */
    void applyEstimated(std::vector<ProteinIdentification>& ids) const;

    /**
    @brief Calculate a linear combination of the area of the difference in estimated vs. empirical (TD) FDR
     and the ROC-N value (AUC up to first N false positives).

    @param ids protein identifications, containing PEP scores annotated with target decoy. Only first run will be evaluated.
    @param pepCutoff up to which PEP should the differences between the two FDRs be calculated
    @param fpCutoff up to which nr. of false positives should the target-decoy AUC be evaluated
    @param diffWeight which weight should the difference get. The ROC-N value gets 1 - this weight.
    */
    double applyEvaluateProteinIDs(const std::vector<ProteinIdentification>& ids, double pepCutoff = 1.0, UInt fpCutoff = 50, double diffWeight = 0.2) const;
    /**
    @brief Calculate a linear combination of the area of the difference in estimated vs. empirical (TD) FDR
     and the ROC-N value (AUC up to first N false positives).

    @param ids protein identifications, containing PEP scores annotated with target decoy.
    @param pepCutoff up to which PEP should the differences between the two FDRs be calculated
    @param fpCutoff up to which nr. of false positives should the target-decoy AUC be evaluated
    @param diffWeight which weight should the difference get. The ROC-N value gets 1 - this weight.
    */
    double applyEvaluateProteinIDs(const ProteinIdentification& ids, double pepCutoff = 1.0, UInt fpCutoff = 50, double diffWeight = 0.2) const;

    /**
    @brief Calculate a linear combination of the area of the difference in estimated vs. empirical (TD) FDR
     and the ROC-N value (AUC up to first N false positives).

    @param score_to_tgt_dec_fraction_pairs extracted scores of protein(group) identifications, containing PEP scores annotated with target decoy fractions. Simple case target=1, decoy=0.
    @param pepCutoff up to which PEP should the differences between the two FDRs be calculated
    @param fpCutoff up to which nr. of false positives should the target-decoy AUC be evaluated
    @param diffWeight which weight should the difference get. The ROC-N value gets 1 - this weight.
    */
    double applyEvaluateProteinIDs(ScoreToTgtDecLabelPairs& score_to_tgt_dec_fraction_pairs, double pepCutoff = 1.0, UInt fpCutoff = 50, double diffWeight = 0.2) const;

    /// simpler reimplementation of the apply function above for PSMs. With charge and identifier info from @p run_info
    void applyBasic(const std::vector<ProteinIdentification> & run_info, std::vector<PeptideIdentification> & ids);

    /// simpler reimplementation of the apply function above for PSMs or peptides.
    void applyBasic(std::vector<PeptideIdentification> & ids, bool higher_score_better, int charge = 0, String identifier = "", bool only_best_per_pep = false);
    /// like applyBasic with "only_best_per_peptide" but it assigns a score to EVERY PSM sharing the peptide sequence with the
    /// best representative. Useful if all hits need to have a peptide score (e.g., for mzTab report). No support for specific charges, runs etc. yet
    void applyBasicPeptideLevel(std::vector<PeptideIdentification> & ids);
    /// like applyBasic with "only_best_per_peptide" but it assigns a score to EVERY PSM sharing the peptide sequence with the
    /// best representative. Useful if all hits need to have a peptide score (e.g., for mzTab report). No support for specific charges, runs etc. yet
    void applyBasicPeptideLevel(ConsensusMap & ids, bool use_unassigned_peptides = true);
    /// simpler reimplementation of the apply function above for peptides in ConsensusMaps.
    void applyBasic(ConsensusMap & cmap, bool use_unassigned_peptides = true);
    /// simpler reimplementation of the apply function above for proteins.
    void applyBasic(ProteinIdentification & id, bool groups_too = true);

    /**
     * @brief  Applies a picked protein FDR.
     * Behaves like a normal target-decoy FDR where only the score of the best protein per
     * target-decoy pair is used. A pair is calculated by checking accession equality after removing the decoy string.
     * If @p decoy_string is empty, we try to guess it. If you set @p decoy_string you should also set @p prefix and
     * say if the string is a prefix (true) or suffix (false).
     * @p groups_too decides if also a (indistinguishable) group-level FDR will be calculated. Here a group score
     * will be taken if not ALL proteins in the group were picked already. Targets preferred.
     */
    void applyPickedProteinFDR(ProteinIdentification& id, String decoy_string = "", bool prefix = true, bool groups_too = true);

    /// calculates the AUC until the first fp_cutoff False positive pep IDs (currently only takes all runs together)
    /// if fp_cutoff = 0, it will calculate the full AUC
    double rocN(const std::vector<PeptideIdentification>& ids, Size fp_cutoff) const;

    /// calculates the AUC until the first fp_cutoff False positive pep IDs (currently only takes all runs together)
    /// if fp_cutoff = 0, it will calculate the full AUC. Restricted to IDs from a specific ID run.
    double rocN(const std::vector<PeptideIdentification>& ids, Size fp_cutoff, const String& identifier) const;

    /// calculates the AUC until the first @p fp_cutoff False positive pep IDs (takes all runs together)
    /// if fp_cutoff = 0, it will calculate the full AUC
    double rocN(const ConsensusMap& ids, Size fp_cutoff, bool include_unassigned_peptides = false) const;

    /// calculates the AUC until the first @p fp_cutoff False positive pep IDs.
    /// if fp_cutoff = 0, it will calculate the full AUC. Restricted to IDs from a specific ID run with @p identifier.
    double rocN(const ConsensusMap& ids, Size fp_cutoff, const String& identifier, bool include_unassigned_peptides = false) const;

    //TODO the next two methods could potentially be merged for speed (they iterate over the same structure)
    //But since they have different cutoff types and it is more generic, I leave it like this.
    /// calculates the area of the difference between estimated and empirical FDR on the fly. Does not store results.
    double diffEstimatedEmpirical(const ScoreToTgtDecLabelPairs& scores_labels, double pepCutoff = 1.0) const;

    /// calculates AUC of empirical FDR up to the first fpCutoff false positives on the fly. Does not store results.
    /// use e.g. fpCutoff = scores_labels.size() for complete AUC
    double rocN(const ScoreToTgtDecLabelPairs& scores_labels, Size fpCutoff = 50) const;

    /**
       @brief Calculate FDR on the level of observation matches (e.g. peptide-spectrum matches) for "general" identification data

       @param id_data Identification data
       @param score_ref Key of the score to use for FDR calculation

       @return Key of the FDR score
    */
    IdentificationData::ScoreTypeRef applyToObservationMatches(IdentificationData& id_data, IdentificationData::ScoreTypeRef score_ref) const;

    /**
     * @brief Finds decoy strings in ProteinIdentification runs
     */
    class DecoyStringHelper
    {
    public:
      /**
       * A result of the findDecoyString function
       */
      struct Result
      {
        bool success; ///< did more than 30% of proteins have the *same* prefix or suffix
        String name; ///< on success, what was the decoy string?
        bool is_prefix; ///< on success, was it a prefix or suffix
      };

      /**
       * @brief Finds the most common decoy string in the accessions of @p proteins. Checks for suffix and prefix and
       * some common decoy strings. Only successful if more than 30% had a common string.
       * @param proteins Input proteins with accessions
       * @return A @struct Result
       */
      static Result findDecoyString(const ProteinIdentification& proteins);
    };
private:

    /// Not implemented
    FalseDiscoveryRate(const FalseDiscoveryRate&);

    /// Not implemented
    FalseDiscoveryRate& operator=(const FalseDiscoveryRate&);

    /// calculates the FDR, given two vectors of scores
    void calculateFDRs_(std::map<double, double>& score_to_fdr, std::vector<double>& target_scores, std::vector<double>& decoy_scores, bool q_value, bool higher_score_better) const;

    /// Helper function for applyToObservationMatches()
    void handleObservationMatch_(
        IdentificationData::ObservationMatchRef match_ref,
        IdentificationData::ScoreTypeRef score_ref,
        std::vector<double>& target_scores,
        std::vector<double>& decoy_scores,
        std::map<IdentificationData::IdentifiedMolecule, bool>& molecule_to_decoy,
        std::map<IdentificationData::ObservationMatchRef, double>& match_to_score) const;

    /// calculates an estimated FDR (based on P(E)Ps) given a vector of score value pairs and fills a map for lookup
    /// in scores_to_FDR
    void calculateEstimatedQVal_(std::map<double, double> &scores_to_FDR,
                                 ScoreToTgtDecLabelPairs &scores_labels,
                                 bool higher_score_better) const;

    /// calculates the FDR with a basic and faster algorithm
    /// Just goes through the sorted scores and counts the number of decoys and targets and annotates the FDR for
    /// this score as it goes. Q-values are optionally annotated by calculating the cumulative minimum in reversed
    /// order afterwards. Since I never understood our other algorithm, I can not explain the difference.
    /// @note Formula used depends on Param "conservative": false -> (D+1)/T, true (e.g. used in Fido) -> (D+1)/(T+D)
    void calculateFDRBasic_(std::map<double,double>& scores_to_FDR, ScoreToTgtDecLabelPairs& scores_labels, bool qvalue, bool higher_score_better) const;

    /// calculates the error area around the x=x line between two consecutive values of expected and actual
    /// i.e. it assumes exp2 > exp1
    double trapezoidal_area_xEqy(double exp1, double exp2, double act1, double act2) const;

    /// calculates the trapezoidal area for a trapezoid with a flat horizontal base e.g. for an AUC
    double trapezoidal_area(double x1, double x2, double y1, double y2) const;
  };

} // namespace OpenMS

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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <boost/unordered_map.hpp>

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
    - true: (D+1)/(T+D) [for comparison with protein level FDR in Fido mostly]
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
    @brief Calculates the FDR of one run from a concatenated sequence DB search

    @param id peptide identifications, containing target and decoy hits
    */
    void apply(std::vector<PeptideIdentification>& id) const;

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

    @param ids protein identifications, containing PEP scores annotated with target decoy. If vector, only first will be evaluated-
    @param pepCutoff up to which PEP should the differences between the two FDRs be calculated
    @param fpCutoff up to which nr. of false positives should the target-decoy AUC be evaluated
    @param diffWeight which weight should the difference get. The ROC-N value gets 1 - this weight.
    */
    double applyEvaluateProteinIDs(const std::vector<ProteinIdentification>& ids, double pepCutoff = 1.0, UInt fpCutoff = 50, double diffWeight = 0.2);
    double applyEvaluateProteinIDs(const ProteinIdentification& ids, double pepCutoff = 1.0, UInt fpCutoff = 50, double diffWeight = 0.2);
    double applyEvaluateProteinIDs(ScoreToTgtDecLabelPairs& score_to_tgt_dec_fraction_pairs, double pepCutoff = 1.0, UInt fpCutoff = 50, double diffWeight = 0.2);

    /// simpler reimplemetation of the apply function above.
    void applyBasic(std::vector<PeptideIdentification> & ids);
    /// simpler reimplemetation of the apply function above for peptides in ConsensusMaps.
    void applyBasic(ConsensusMap & cmap, bool use_unassigned_peptides = true);
    /// simpler reimplemetation of the apply function above for proteins.
    void applyBasic(ProteinIdentification & id, bool groups_too = true);

    /// calculates the AUC until the first fp_cutoff False positive pep IDs (currently only takes all runs together)
    /// if fp_cutoff = 0, it will calculate the full AUC
    double rocN(const std::vector<PeptideIdentification>& ids, Size fp_cutoff) const;

    /// calculates the AUC until the first fp_cutoff False positive pep IDs (currently only takes all runs together)
    /// if fp_cutoff = 0, it will calculate the full AUC. Restricted to IDs from a specific ID run.
    double rocN(const std::vector<PeptideIdentification>& ids, Size fp_cutoff, const String& identifier) const;

    /// calculates the AUC until the first fp_cutoff False positive pep IDs (currently only takes all runs together)
    /// if fp_cutoff = 0, it will calculate the full AUC
    double rocN(const ConsensusMap& ids, Size fp_cutoff) const;

    /// calculates the AUC until the first fp_cutoff False positive pep IDs (currently only takes all runs together)
    /// if fp_cutoff = 0, it will calculate the full AUC. Restricted to IDs from a specific ID run.
    double rocN(const ConsensusMap& ids, Size fp_cutoff, const String& identifier) const;

    //TODO the next two methods could potentially be merged for speed (they iterate over the same structure)
    //But since they have different cutoff types and it is more generic, I leave it like this.
    /// calculates the area of the difference between estimated and empirical FDR on the fly. Does not store results.
    double diffEstimatedEmpirical(const ScoreToTgtDecLabelPairs& scores_labels, double pepCutoff = 1.0) const;

    /// calculates AUC of empirical FDR up to the first fpCutoff false positives on the fly. Does not store results.
    /// use e.g. fpCutoff = scores_labels.size() for complete AUC
    double rocN(const ScoreToTgtDecLabelPairs& scores_labels, Size fpCutoff = 50) const;

    /**
       @brief Calculate FDR on the level of molecule-query matches (e.g. peptide-spectrum matches) for "general" identification data

       @param id_data Identification data
       @param score_key Key of the score to use for FDR calculation

       @return Key of the FDR score
    */
    IdentificationData::ScoreTypeRef applyToQueryMatches(IdentificationData& id_data, IdentificationData::ScoreTypeRef score_ref) const;


private:

    /// Not implemented
    FalseDiscoveryRate(const FalseDiscoveryRate&);

    /// Not implemented
    FalseDiscoveryRate& operator=(const FalseDiscoveryRate&);

    /// calculates the FDR, given two vectors of scores
    void calculateFDRs_(std::map<double, double>& score_to_fdr, std::vector<double>& target_scores, std::vector<double>& decoy_scores, bool q_value, bool higher_score_better) const;

    /// Helper function for applyToQueryMatches()
    void handleQueryMatch_(
        IdentificationData::QueryMatchRef match_ref,
        IdentificationData::ScoreTypeRef score_ref,
        std::vector<double>& target_scores,
        std::vector<double>& decoy_scores,
        std::map<IdentificationData::IdentifiedMoleculeRef, bool>& molecule_to_decoy,
        std::map<IdentificationData::QueryMatchRef, double>& match_to_score) const;

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


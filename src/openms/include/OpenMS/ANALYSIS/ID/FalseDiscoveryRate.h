// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
    "target_decoy" metavalue e.g. annotated by PeptideIndexer

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

    /// simpler reimplemetation of the apply function above.
    void applyBasic(std::vector<PeptideIdentification> & ids);
    /// @todo groups and proteins not supported on CMaps yet. WIP
    void applyBasic(ConsensusMap & cmap, bool groups, bool proteins, bool peptides);
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

    template<typename T> struct IsIDType
    {
      static bool const value = std::is_same<T, PeptideIdentification>::value || std::is_same<T, ProteinIdentification>::value;
    };

    template<typename T> struct IsHitType
    {
      static bool const value = std::is_same<T, PeptideHit>::value || std::is_same<T, ProteinHit>::value;
    };

    /// Used to collect data from the ID structures with the original score as first and
    /// target decoy annotation as second member of the pair. Target = true.
    /// Target+decoy for peptides = target. Protein groups with at least one target = target.
    typedef std::vector<std::pair<double,bool>> ScoreToTgtDecLabelPairs;

    /**
     * \defgroup getScoresFunctions Get scores from ID structures for FDR
     * @brief  Fills the scores_labels vector from an ID data structure
     * @param  scores_labels Pairs of scores and boolean target decoy labels to be filled. target = true.
     *
     * Just use the one you need.
     * @{
     */


    //TODO could be done with set of target accessions, too
    //TODO even better: store nr targets and nr decoys when creating the groups!
    //TODO alternative scoring is possible, too (e.g. ratio of tgts vs decoys),
    // this requires templatization of the scores_labels vector though
    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const std::vector<ProteinIdentification::ProteinGroup>& grps,
        const std::unordered_set<std::string>& decoy_accs) const;

    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const ProteinIdentification & id) const
    {

      scores_labels.reserve(scores_labels.size() + id.getHits().size());
      std::transform(id.getHits().begin(), id.getHits().end(),
                     std::back_inserter(scores_labels),
                     [](const ProteinHit& hit)
                     {
                       checkTDAnnotation_(hit);
                       return std::make_pair<double,bool>(hit.getScore(), getTDLabel_(hit));
                     }
                    );

    }

    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const PeptideIdentification & id, bool all_hits, int charge, const String& identifier) const
    {
      if (id.getIdentifier() == identifier)
      {
        getScores_(scores_labels, id, all_hits, charge);
      }
    }


    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const PeptideIdentification & id, bool all_hits, const String& identifier) const
    {
      if (id.getIdentifier() == identifier)
      {
        getScores_(scores_labels, id, all_hits);
      }
    }

    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const PeptideIdentification & id, int charge, const String& identifier) const
    {
      if (id.getIdentifier() == identifier)
      {
        getScores_(scores_labels, id, charge);
      }
    }

    template<typename IDType, typename std::enable_if<IsIDType<IDType>::value>::type* = nullptr>
    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const IDType & id, const String& identifier) const
    {
      if (id.getIdentifier() == identifier)
      {
        getScores_(scores_labels, id);
      }
    }

    template<class ...Args>
    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const PeptideIdentification & id,
        bool all_hits,
        Args&& ... args) const
    {
      if (all_hits)
      {
        for (const PeptideHit &hit : id.getHits())
        {
          getScores_(scores_labels, hit, std::forward<Args>(args)...);
        }
      }
      else
      {
        //TODO for speed I assume that they are sorted and first = best.
        //id.sort();
        const PeptideHit& hit = id.getHits()[0];
        getScores_(scores_labels, hit, std::forward<Args>(args)...);
      }
    }

    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const PeptideHit & hit,
        int charge) const
    {
        if (charge == hit.getCharge())
        {
          checkTDAnnotation_(hit);
          scores_labels.emplace_back(hit.getScore(), getTDLabel_(hit));
        }
    }

    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const PeptideIdentification & id,
        int charge) const
    {
      for (const PeptideHit &hit : id.getHits())
      {
        getScores_(scores_labels, hit, charge);
      }
    }

    template<typename HitType, typename std::enable_if<IsHitType<HitType>::value>::type* = nullptr>
    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const HitType & hit) const
    {
        checkTDAnnotation_(hit);
        scores_labels.emplace_back(hit.getScore(), getTDLabel_(hit));
    }

    template<typename IDType, typename std::enable_if<IsIDType<IDType>::value>::type* = nullptr>
    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const IDType & id) const
    {
      for (const typename IDType::HitType &hit : id.getHits())
      {
        getScores_(scores_labels, hit);
      }
    }

    template<class ...Args>
    void getScores_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const std::vector<PeptideIdentification> & ids,
        Args&& ... args) const
    {
      for (const auto& id : ids)
      {
        getScores_(scores_labels, id, std::forward<Args>(args)...);
      }
    }
    /** @} */

    /**
     * @brief Helper for getting scores in ConsensusMaps
     * @todo allow FeatureMap?
     */
    template<class ...Args>
    void getPeptideScoresFromMap_(
        ScoreToTgtDecLabelPairs& scores_labels,
        const ConsensusMap & cmap, Args&& ... args) const
    {
      std::function<void (const PeptideIdentification &)> f =
          [&, this](const PeptideIdentification& id) -> void {this->getScores_(scores_labels, id, std::forward<Args>(args)...);};
      cmap.applyFunctionOnPeptideIDs(f);
    }

    /**
     * @brief For peptide hits, a hit is considered target also if it maps to both
     * a target and a decoy protein (i.e. "target+decoy") as value in the
     * "target_decoy" metavalue e.g. annotated by PeptideIndexer
     */
    static bool getTDLabel_(const MetaInfoInterface& idOrHit)
    {
      return std::string(idOrHit.getMetaValue("target_decoy"))[0] == 't';
    }

    /**
     * \defgroup setScoresFunctions Sets FDRs/qVals
     * @brief  Sets FDRs/qVals from a scores_to_FDR map in the ID data structures
     * @param  scores_to_FDR Maps original score to calculated FDR or q-Value
     * @param  score_type FDR or q-Value
     * @param  higher_better should usually be false @todo remove?
     *
     * Just use the one you need.
     * @{
     */

    template <typename IDType, class ...Args>
    void setScores_(const std::map<double,double>& scores_to_FDR, std::vector<IDType> & ids, const std::string& score_type,
                    bool higher_better, Args& ... args) const
    {
      for (auto& id : ids)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better, &args...);
      }
    }

    template <typename IDType>
    String setScoreType_(IDType& id, const std::string& score_type,
                    bool higher_better) const
    {
      String old_score_type = id.getScoreType() + "_score";
      id.setScoreType(score_type);
      id.setHigherScoreBetter(higher_better);
      return old_score_type;
    }

    template <typename IDType>
    void setScores_(const std::map<double,double>& scores_to_FDR, IDType & id, const std::string& score_type,
        bool higher_better, bool keep_decoy) const
    {
      String old_score_type = setScoreType_(id, score_type, higher_better);

      if (keep_decoy) //in-place set scores
      {
        setScores_(scores_to_FDR, id, old_score_type);
      }
      else
      {
        setScoresAndRemoveDecoys_(scores_to_FDR, id, old_score_type);
      }
    }

    template <typename IDType>
    void setScores_(const std::map<double,double>& scores_to_FDR, IDType & id,
                    const String& old_score_type) const
    {
      std::vector<typename IDType::HitType>& hits = id.getHits();
      for (auto &hit : hits)
      {
        setScore_(scores_to_FDR, hit, old_score_type);
      }
    }

    template <typename IDType, class ...Args>
    void setScoresAndRemoveDecoys_(const std::map<double,double>& scores_to_FDR, IDType & id,
                    const String& old_score_type, Args ... args) const
    {
      std::vector<typename IDType::HitType>& hits = id.getHits();
      std::vector<typename IDType::HitType> new_hits;
      new_hits.reserve(hits.size());
      for (auto &hit : hits)
      {
        setScoreAndMoveIfTarget_( scores_to_FDR, hit, old_score_type, new_hits, args...);
      }
      hits.swap(new_hits);
    }

    template <typename HitType>
    void setScore_(const std::map<double,double>& scores_to_FDR, HitType & hit, const std::string& old_score_type) const
    {
      hit.setMetaValue(old_score_type, hit.getScore());
      hit.setScore(scores_to_FDR.lower_bound(hit.getScore())->second);
    }

    template <typename IDType>
    void setScores_(const std::map<double,double>& scores_to_FDR, IDType & id, const std::string& score_type,
                    bool higher_better) const
    {
      String old_score_type = setScoreType_(id, score_type, higher_better);
      setScores_(scores_to_FDR, id, old_score_type);
    }

    void setScores_(const std::map<double,double>& scores_to_FDR, PeptideIdentification& id, const std::string& score_type,
                    bool higher_better, bool keep_decoy, int charge) const
    {
      String old_score_type = setScoreType_(id, score_type, higher_better);
      if (keep_decoy) //in-place set scores
      {
        setScores_(scores_to_FDR, id, old_score_type, charge);
      }
      else
      {
        setScoresAndRemoveDecoys_<PeptideIdentification>(scores_to_FDR, id, old_score_type, charge);
      }
    }

    void setScores_(const std::map<double,double>& scores_to_FDR, PeptideIdentification & id, const std::string& score_type,
                    bool higher_better, bool keep_decoy, int charge, const String& identifier) const
    {
      if (id.getIdentifier() == identifier)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better, keep_decoy, charge);
      }
    }

    template <typename IDType>
    void setScores_(const std::map<double,double>& scores_to_FDR, IDType & id, const std::string& score_type,
                    bool higher_better, bool keep_decoy, const String& identifier) const
    {
      if (id.getIdentifier() == identifier)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better, keep_decoy);
      }
    }

    void setScores_(const std::map<double,double>& scores_to_FDR, PeptideIdentification & id, const std::string& score_type,
                    bool higher_better, int charge, const String& identifier) const
    {
      if (id.getIdentifier() == identifier)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better, charge);
      }
    }

    template <typename IDType>
    void setScores_(const std::map<double,double>& scores_to_FDR, IDType & id, const std::string& score_type,
                    bool higher_better, const String& identifier) const
    {
      if (id.getIdentifier() == identifier)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better);
      }
    }

    //TODO could also get a keep_decoy flag when we define what a "decoy group" is -> keep all always for now
    void setScores_(
        const std::map<double,double>& scores_to_FDR,
        std::vector<ProteinIdentification::ProteinGroup>& grps,
        const std::string& score_type,
        bool higher_better) const;

    /** @} */

    /**
     * @brief Used when keep_decoy_peptides or proteins is false
     * @tparam HitType Protein or PeptideHit
     * @param scores_to_FDR map from original score to FDR/qVal
     * @param hit the Hit itself
     * @param old_score_type to save it in metavalue
     * @param new_hits where to move if target (i.e. target or target+decoy)
     */
    template <typename HitType>
    void setScoreAndMoveIfTarget_(const std::map<double,double>& scores_to_FDR, HitType & hit, const std::string& old_score_type,
                                  std::vector<HitType>& new_hits) const
    {
      const String& target_decoy(hit.getMetaValue("target_decoy"));
      if (target_decoy[0] == 't')
      {
        hit.setMetaValue(old_score_type, hit.getScore());
        hit.setScore(scores_to_FDR.at(hit.getScore()));
        new_hits.push_back(std::move(hit));
      } // else do not move over
    }

     /**
     * @brief Used when keep_decoy_peptides is false and charge states are considered
     * @param scores_to_FDR map from original score to FDR/qVal
     * @param hit the PeptideHit itself
     * @param old_score_type to save it in metavalue
     * @param new_hits where to move if target (i.e. target or target+decoy)
     * @param charge If only peptides with charge X are currently considered
     */
    void setScoreAndMoveIfTarget_(const std::map<double,double>& scores_to_FDR, PeptideHit & hit, const std::string& old_score_type,
                                  std::vector<PeptideHit>& new_hits, int charge) const
    {
      if (charge == hit.getCharge())
      {
        const String &target_decoy(hit.getMetaValue("target_decoy"));
        if (target_decoy[0] == 't')
        {
          hit.setMetaValue(old_score_type, hit.getScore());
          hit.setScore(scores_to_FDR.at(hit.getScore()));
          new_hits.push_back(std::move(hit));
        } // else do not move over
      }
      else // different charge, move over unchanged to process later at correct charge.
      {
        new_hits.push_back(std::move(hit));
      }
    }

    /**
     * @brief Helper for applying set Scores on ConsensusMaps
     * @tparam Args optional additional arguments (charge, run ID)
     * @param scores_to_FDR maps original scores to FDR
     * @param cmap the ConsensusMap
     * @param score_type FDR or q-Value
     * @param higher_better usually false
     * @param keep_decoy read from Param object
     * @param args optional additional arguments (int charge, string run ID)
     */
    template <class ...Args>
    void setPeptideScoresForMap_(const std::map<double,double>& scores_to_FDR, ConsensusMap& cmap,
                                 const std::string& score_type, bool higher_better, bool keep_decoy, Args&& ... args) const
    {
      std::function<void (PeptideIdentification &)> f =
          [&,this](PeptideIdentification& id) -> void {this->setScores_(scores_to_FDR, id, score_type, higher_better, keep_decoy, std::forward<Args>(args)...);};
      cmap.applyFunctionOnPeptideIDs(f);
    }

    /**
     * @brief To check the metavalues before we do anything
     * @param idOrHit
     * @throws Exception::MissingInformation if it does not exist
     */
    static void checkTDAnnotation_(const MetaInfoInterface & idOrHit)
    {
        if (!idOrHit.metaValueExists("target_decoy"))
        {
          throw Exception::MissingInformation(__FILE__,
                                              __LINE__,
                                              OPENMS_PRETTY_FUNCTION,
                                              "Meta value 'target_decoy' does not exist in all ProteinHits! Reindex the idXML file with 'PeptideIndexer'");
        }
    }

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
    void calculateFDRBasic_(std::map<double,double>& scores_to_FDR, ScoreToTgtDecLabelPairs& scores_labels, bool qvalue, bool higher_score_better);

    //TODO the next two methods could potentially be merged for speed (they iterate over the same structure)
    //But since they have different cutoff types and it is more generic, I leave it like this.
    /// calculates the area of the difference between estimated and empirical FDR on the fly. Does not store results.
    double diffEstimatedEmpirical_(const ScoreToTgtDecLabelPairs& scores_labels, double pepCutoff = 1.0);

    /// calculates AUC of empirical FDR up to the first fpCutoff false positives on the fly. Does not store results.
    /// use e.g. fpCutoff = scores_labels.size() for complete AUC
    double rocN_(ScoreToTgtDecLabelPairs const &scores_labels, Size fpCutoff = 50) const;

    /// calculates the error area around the x=x line between two consecutive values of expected and actual
    /// i.e. it assumes exp2 > exp1
    double trapezoidal_area_xEqy(double exp1, double exp2, double act1, double act2) const;

    /// calculates the trapezoidal area for a trapezoid with a flat horizontal base e.g. for an AUC
    double trapezoidal_area(double x1, double x2, double y1, double y2) const;

  };

} // namespace OpenMS


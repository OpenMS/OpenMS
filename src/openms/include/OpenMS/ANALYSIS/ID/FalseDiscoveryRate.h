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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Calculates an FDR from identifications

        Either two runs of forward and decoy database identification or
        one run containing both (with marks) can be used to annotate
        each of the peptide hits with a FDR.

    Also q-values can be reported instead of p-values.
    q-values are basically only adjusted p-values, also ranging from 0 to 1, with lower values being preferable.
    When looking at the list of hits ordered by q-values, then a hit with q-value of @em x means that there is an
    @em x*100 percent chance that all hits with a q-value <= @em x are a false positive hit.

        @todo implement combined searches properly (Andreas)
        @improvement implement charge state separated fdr/q-values (Andreas)

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
    void apply(std::vector<PeptideIdentification> & fwd_ids, std::vector<PeptideIdentification> & rev_ids) const;

    /**
    @brief Calculates the FDR of one run from a concatenated sequence db search

    @param id peptide identifications, containing target and decoy hits
    */
    void apply(std::vector<PeptideIdentification> & id) const;

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
    @brief Calculate the FDR based on PEPs pr PPs (if present) and modifies the IDs inplace

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

    void applyBasic(std::vector<PeptideIdentification> & ids);
    void applyBasic(ProteinIdentification & id);


private:

    ///Not implemented
    FalseDiscoveryRate(const FalseDiscoveryRate &);

    ///Not implemented
    FalseDiscoveryRate & operator=(const FalseDiscoveryRate &);

    //TODO we could add identifier here. If we need to combine runs.
    void getScores_(
      std::vector<std::pair<double,bool>>& scores_labels, 
      const ProteinIdentification & id) const;

    void getScores_(
      std::vector<std::pair<double,bool>>& scores_labels, 
      const std::vector<PeptideIdentification> & ids, 
      bool all_hits, 
      int charge, String identifier) const;

    void getScores_(
      std::vector<std::pair<double,bool>>& scores_labels, 
      const std::vector<PeptideIdentification> & targets, 
      const std::vector<PeptideIdentification> & decoys, 
      bool all_hits, 
      int charge, 
      const String& identifier) const;

    void setScores_(
      const std::map<double,double>& scores_to_FDR, 
      std::vector<PeptideIdentification> & id, 
      const std::string& score_type, 
      bool higher_better) const;

    template <typename IDType>
    void setScores_(const std::map<double,double>& scores_to_FDR, IDType & id, const std::string& score_type, bool higher_better) const
    {
      String old_score_type = id.getScoreType() + "_score";
      id.setScoreType(score_type);
      id.setHigherScoreBetter(higher_better);
      for (auto& hit : id.getHits())
      {
        double old_score = hit.getScore();
        hit.setScore(scores_to_FDR.lower_bound(hit.getScore())->second);
        hit.setMetaValue(old_score_type, old_score);
      }
    }

    template <typename IDType>
    void checkTDAnnotation_ (const IDType & id) const
    {
      for (auto const& hit : id.getHits())
      {
        if (!hit.metaValueExists("target_decoy"))
        {
          throw Exception::MissingInformation(__FILE__,
                                              __LINE__,
                                              OPENMS_PRETTY_FUNCTION,
                                              "Meta value 'target_decoy' does not exist in all ProteinHits! Reindex the idXML file with 'PeptideIndexer'");
        }
      }
    }

    template <typename HitType>
    struct GetLabelFunctor: std::function<bool(const HitType&)>
    {
      bool operator() (const HitType& hit)
      {
          if (!hit.metaValueExists("target_decoy"))
          {
            throw Exception::MissingInformation(__FILE__,
                                                __LINE__,
                                                OPENMS_PRETTY_FUNCTION,
                                                "Meta value 'target_decoy' does not exist in all ProteinHits! Reindex the idXML file with 'PeptideIndexer'");
          }
          else
          {
            return std::string(hit.getMetaValue("target_decoy"))[0] == 't';
          }
      }
    };

    template <typename HitType>
    struct TrueFunctor: std::function<bool(const HitType&)>
    {
      bool operator() (const HitType& /*hit*/)
      {
        return true;
      }
    };

    template <typename HitType>
    struct FalseFunctor: std::function<bool(const HitType&)>
    {
      bool operator() (const HitType& /*hit*/)
      {
        return false;
      }
    };


    template <typename HitType>
    std::pair<double,bool> getScoreLabel_(const HitType& hit, std::function<bool(const HitType&)> fun) const
    {
      return std::make_pair(hit.getScore(), fun(hit));
    }


    /// calculates the fdr given two vectors of scores and fills a map for lookup in scores_to_FDR
    void calculateFDRs_(Map<double, double>& score_to_fdr, std::vector<double>& target_scores, std::vector<double>& decoy_scores, bool q_value, bool higher_score_better) const;

    /// calculates an estimated FDR (based on P(E)Ps) given a vector of score value pairs and fills a map for lookup
    /// in scores_to_FDR
    void calculateEstimatedQVal_(std::map<double, double> &scores_to_FDR,
                                 std::vector<std::pair<double, bool>> &scores_labels,
                                 bool higher_score_better) const;

    /// calculates the FDR with a basic and faster algorithm
    void calculateFDRBasic_(std::map<double,double>& scores_to_FDR, std::vector<std::pair<double,bool>>& scores_labels, bool qvalue, bool higher_score_better);

    //TODO the next two methods could potentially be merged for speed (they iterate over the same structure)
    //But since they have different cutoff types and it is more generic, I leave it like this.
    /// calculates the area of the difference between estimated and  empirical FDR on the fly. Does not store results.
    double diffEstimatedEmpirical_(const std::vector<std::pair<double, bool>>& scores_labels, double pepCutoff = 1.0);
    /// calculates AUC of empirical FDR up to the first fpCutoff false positives on the fly. Does not store results.
    double rocN_(std::vector<std::pair<double, bool>> const &scores_labels, UInt fpCutoff = 50);

    /// calculates the error area around the x=x line between two consecutive values of expected and actual
    /// i.e. it assumes exp2 > exp1
    double trapezoidal_area_xEqy(double exp1, double exp2, double act1, double act2);

    /// calculates the trapezoidal area for a trapezoid with a flat horizontal base e.g. for an AUC
    double trapezoidal_area(double x1, double x2, double y1, double y2);


  };

} // namespace OpenMS


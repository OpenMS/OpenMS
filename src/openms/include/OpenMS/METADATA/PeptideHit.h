// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PEPTIDEHIT_H
#define OPENMS_METADATA_PEPTIDEHIT_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideEvidence.h>

namespace OpenMS
{
  /**
    @brief Representation of a peptide hit

    It contains the fields score, score_type, rank, and sequence.

        @ingroup Metadata
  */
  class OPENMS_DLLAPI PeptideHit :
    public MetaInfoInterface
  {
public:

    /// @name Comparators for PeptideHit and ProteinHit
    //@{
    /// Greater predicate for scores of hits
    class OPENMS_DLLAPI ScoreMore
    {
public:
      template <typename Arg>
      bool operator()(const Arg& a, const Arg& b)
      {
        return a.getScore() > b.getScore();
      }

    };

    /// Lesser predicate for scores of hits
    class OPENMS_DLLAPI ScoreLess
    {
public:
      template <typename Arg>
      bool operator()(const Arg& a, const Arg& b)
      {
        return a.getScore() < b.getScore();
      }

    };

    /// Lesser predicate for scores of hits
    class OPENMS_DLLAPI RankLess
    {
public:
      template <typename Arg>
      bool operator()(const Arg& a, const Arg& b)
      {
        return a.getRank() < b.getRank();
      }

    };
    //@}

    /// Analysis Result (containing search engine / prophet results)
    class OPENMS_DLLAPI PepXMLAnalysisResult
    {
public:
      String score_type; // e.g. peptideprophet / interprophet
      bool higher_is_better; // is higher score better ?
      double main_score; // posterior probability for example
      std::map<String, double> sub_scores; /// additional scores attached to the original, aggregated score

      bool operator==(const PepXMLAnalysisResult& rhs) const
      {
        return score_type == rhs.score_type 
          && higher_is_better == rhs.higher_is_better
          && main_score == rhs.main_score
          && sub_scores == rhs.sub_scores;
      }

      PepXMLAnalysisResult& operator=(const PepXMLAnalysisResult& source)
      {
        if (this == &source) return *this;
        score_type = source.score_type;
        higher_is_better = source.higher_is_better;
        main_score = source.main_score;
        sub_scores = source.sub_scores;
        return *this;
      }

    };

    /** @name Constructors and Destructor */
    //@{
    /// default constructor
    PeptideHit();

    /// values constructor
    PeptideHit(double score,
               UInt rank,
               Int charge,
               const AASequence& sequence);

    /// copy constructor
    PeptideHit(const PeptideHit& source);

    /// destructor
    virtual ~PeptideHit();
    //@}

    /// assignment operator
    PeptideHit& operator=(const PeptideHit& source);

    /// Equality operator
    bool operator==(const PeptideHit& rhs) const;

    /// Inequality operator
    bool operator!=(const PeptideHit& rhs) const;

    /**	@name Accessors
    */
    //@{
    /// returns the peptide sequence without trailing or following spaces
    const AASequence& getSequence() const;

    /// sets the peptide sequence
    void setSequence(const AASequence& sequence);

    /// returns the charge of the peptide
    Int getCharge() const;

    /// sets the charge of the peptide
    void setCharge(Int charge);

    /// returns information on peptides (potentially) identified by this PSM
    const std::vector<PeptideEvidence>& getPeptideEvidences() const;

    /// set information on peptides (potentially) identified by this PSM
    void setPeptideEvidences(const std::vector<PeptideEvidence>& peptide_evidences);

    /// adds information on a peptide that is (potentially) identified by this PSM
    void addPeptideEvidence(const PeptideEvidence& peptide_evidence);

    /// returns the PSM score
    double getScore() const;

    /// sets the PSM score
    void setScore(double score);

    /// set information on (search engine) sub scores associated with this PSM
    void setAnalysisResults(std::vector<PepXMLAnalysisResult> aresult);

    /// add information on (search engine) sub scores associated with this PSM
    void addAnalysisResults(PepXMLAnalysisResult aresult);

    /// returns information on (search engine) sub scores associated with this PSM
    const std::vector<PepXMLAnalysisResult>& getAnalysisResults() const;

    /// returns the PSM rank
    UInt getRank() const;

    /// sets the PSM rank
    void setRank(UInt newrank);
    //@}

    /// extracts the set of non-empty protein accessions from peptide evidences
    std::set<String> extractProteinAccessions() const;

protected:
    AASequence sequence_;

    /// the score of the peptide hit
    double score_;

    /// additional scores attached to the original, aggregated score
    std::vector<PepXMLAnalysisResult>* analysis_results_;

    /// the position(rank) where the hit appeared in the hit list
    UInt rank_;

    /// the charge of the peptide
    Int charge_;

    /// information on the potential peptides observed through this PSM.
    std::vector<PeptideEvidence> peptide_evidences_;
  };

} // namespace OpenMS

#endif // OPENMS_METADATA_PEPTIDEHIT_H

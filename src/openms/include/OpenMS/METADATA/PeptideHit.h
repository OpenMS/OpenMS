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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

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

    /**
   * @brief Contains annotations of a peak

      The mz and intensity values contain the same information as a spectrum
      would have about the peaks, and can be used to map the additional
      information to the correct peak or reconstruct the annotated spectrum.
      Additionally the charge of the peak and an arbitrary string annotation
      can be stored.

      The string annotation can be e.g. a fragment type like "y3".
      This information can be used e.g. to label peaks in TOPPView.

      The specific application in OpenProXL uses a more complex syntax to
      define the larger number of different ion types found in XL-MS data.

      In the example "[alpha|ci$y3-H2O-NH3]" "alpha" or "beta" determines on
      which of the two peptides the fragmentation occurred, "ci" or "xi"
      determines whether the cross-link and with it the other peptide is
      contained in the fragment, and the last part is the ion type with the
      fragmentation position (index) and losses.  The separators "|" and "$"
      are used to separate the parts easily when parsing the annotation.

   */
  struct PeakAnnotation
  {
    String annotation = "";  // e.g. [alpha|ci$y3-H2O-NH3]
    int charge = 0;
    double mz = -1.;
    double intensity = 0.;

    bool operator<(const PeptideHit::PeakAnnotation& other) const
    {
      // sensible to sort first by m/z and charge
      if (mz < other.mz)
      {
        return true;
      }
      else if (mz > other.mz)
      {
        return false;
      }

      if (charge < other.charge)
      {
        return true;
      }
      else if (charge > other.charge)
      {
        return false;
      }

      if (annotation < other.annotation)
      {
        return true;
      }
      else if (annotation > other.annotation)
      {
        return false;
      }

      if (intensity < other.intensity)
      {
        return true;
      }
      else if (intensity > other.intensity)
      {
        return false;
      }

      return false;
    }

    bool operator==(const PeptideHit::PeakAnnotation& other) const
    {
      if (charge != other.charge || mz != other.mz ||
          intensity != other.intensity || annotation != other.annotation) return false;
      return true;
    }

    static void writePeakAnnotationsString_(String& annotation_string, std::vector<PeptideHit::PeakAnnotation> annotations)
    {
      if (annotations.empty()) { return; }

      // sort by mz, charge, ...
      stable_sort(annotations.begin(), annotations.end());

      String val;
      for (auto& a : annotations)
      {
        annotation_string += String(a.mz) + "," + String(a.intensity) + "," + String(a.charge) + "," + String(a.annotation).quote();
        if (&a != &annotations.back()) { annotation_string += "|"; }
      }
    }

  };

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


    /// Lesser predicate for (modified) sequence of hits
    class OPENMS_DLLAPI SequenceLessComparator
    {
      template <typename Arg>
      bool operator()(const Arg& a, const Arg& b)
      {
        if (a.getSequence().toString() < b.getSequence().toString()) return true;
        return false;
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
    };

    /** @name Constructors and Assignment
    */
    //@{
    /// Default constructor
    PeptideHit();
    /// Values constructor that copies sequence
    PeptideHit(double score,
               UInt rank,
               Int charge,
               const AASequence& sequence);
    /// Values constructor that moves sequence R-value
    PeptideHit(double score,
               UInt rank,
               Int charge,
               AASequence&& sequence);
    /// Copy constructor
    PeptideHit(const PeptideHit& source);
    /// Move constructor
    PeptideHit(PeptideHit&&) noexcept;
    /// Destructor
    virtual ~PeptideHit();

    /// Assignment operator
    PeptideHit& operator=(const PeptideHit& source);
    /// Move assignment operator
    PeptideHit& operator=(PeptideHit&&) noexcept;
    //@}

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

    /// sets the peptide sequence
    void setSequence(AASequence&& sequence);

    /// returns the charge of the peptide
    Int getCharge() const;

    /// sets the charge of the peptide
    void setCharge(Int charge);

    /// returns information on peptides (potentially) identified by this PSM
    const std::vector<PeptideEvidence>& getPeptideEvidences() const;

    /// set information on peptides (potentially) identified by this PSM
    void setPeptideEvidences(const std::vector<PeptideEvidence>& peptide_evidences);

    void setPeptideEvidences(std::vector<PeptideEvidence>&& peptide_evidences);

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

    /// returns the fragment annotations
    std::vector<PeptideHit::PeakAnnotation> getPeakAnnotations() const;

    /// sets the fragment annotations
    void setPeakAnnotations(std::vector<PeptideHit::PeakAnnotation> frag_annotations);

    //@}

    /// extracts the set of non-empty protein accessions from peptide evidences
    std::set<String> extractProteinAccessionsSet() const;

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

    /// annotations of fragments in the corresponding spectrum
    std::vector<PeptideHit::PeakAnnotation> fragment_annotations_;
  };

} // namespace OpenMS

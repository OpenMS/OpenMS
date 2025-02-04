// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <iosfwd>
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
      return std::tie(mz, charge, annotation, intensity) < std::tie(other.mz, other.charge, other.annotation, other.intensity);
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
      String score_type; /// e.g. peptideprophet / interprophet
      bool higher_is_better{}; /// is higher score better ?
      double main_score{}; /// posterior probability for example
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
    /// returns the peptide sequence
    const AASequence& getSequence() const;

    /// returns the mutable peptide sequence
    AASequence& getSequence();

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
    void addAnalysisResults(const PepXMLAnalysisResult& aresult);

    /// returns information on (search engine) sub scores associated with this PSM
    const std::vector<PepXMLAnalysisResult>& getAnalysisResults() const;

    /// returns the PSM rank
    UInt getRank() const;

    /// sets the PSM rank
    void setRank(UInt newrank);

    /// returns the fragment annotations
    std::vector<PeptideHit::PeakAnnotation>& getPeakAnnotations();
    const std::vector<PeptideHit::PeakAnnotation>& getPeakAnnotations() const;


    /// sets the fragment annotations
    void setPeakAnnotations(std::vector<PeptideHit::PeakAnnotation> frag_annotations);

    //@}

    /// extracts the set of non-empty protein accessions from peptide evidences
    std::set<String> extractProteinAccessionsSet() const;

protected:
    AASequence sequence_;

    /// the score of the peptide hit
    double score_{};

    /// additional scores attached to the original, aggregated score
    std::vector<PepXMLAnalysisResult>* analysis_results_;

    /// the position(rank) where the hit appeared in the hit list
    UInt rank_{};

    /// the charge of the peptide
    Int charge_{};

    /// information on the potential peptides observed through this PSM.
    std::vector<PeptideEvidence> peptide_evidences_;

    /// annotations of fragments in the corresponding spectrum
    std::vector<PeptideHit::PeakAnnotation> fragment_annotations_;
  };

  /// Stream operator
  OPENMS_DLLAPI std::ostream& operator<< (std::ostream& stream, const PeptideHit& hit);
} // namespace OpenMS

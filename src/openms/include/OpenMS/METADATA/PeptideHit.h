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
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideEvidence.h>

namespace OpenMS
{
  class PeptideHit;
  using SpectrumMatch = PeptideHit; // better name that might become the default in future version

  /**
    @brief Represents a single spectrum match (candidate) for a specific tandem mass spectrum (MS/MS).

    Stores the primary information about a potential match, including:
    - The sequence (potentially with modifications) using AASequence.
    - The primary score assigned by the identification algorithm (e.g., search engine).
    - The rank of this hit compared to other hits for the same spectrum (stored as a meta value with key "rank").
    - The precursor charge state assumed for this match.
    - Evidence linking the peptide sequence to specific protein sequences (PeptideEvidence).
    - Optional annotations mapping fragment ions in the MS/MS spectrum to interpretations (PeakAnnotation).
    - Optional secondary scores from post-processing tools (PepXMLAnalysisResult).

    Objects are typically contained within a PeptideIdentification object, which represents
    all hits found for a single spectrum. Inherits from MetaInfoInterface, allowing
    arbitrary metadata (key-value pairs) to be attached.

    @deprecated Use SpectrumMatch instead. PeptideHit may be removed in a future OpenMS version.

    @see PeptideIdentification, AASequence, PeptideEvidence, PeakAnnotation, PepXMLAnalysisResult, MetaInfoInterface

    @ingroup Metadata
  */
  class OPENMS_DLLAPI PeptideHit :
    public MetaInfoInterface
  {
public:
    /// Enum for target/decoy annotation
    enum class TargetDecoyType
    {
      TARGET,       ///< Only matches target proteins
      DECOY,        ///< Only matches decoy proteins
      TARGET_DECOY, ///< Matches BOTH target and decoy proteins
      UNKNOWN       ///< Target/decoy status is unknown (meta value not set)
    };

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
  struct OPENMS_DLLAPI PeakAnnotation
  {
    String annotation = "";  // e.g. [alpha|ci$y3-H2O-NH3]
    int charge = 0;
    double mz = -1.;
    double intensity = 0.;

    bool operator<(const PeptideHit::PeakAnnotation& other) const;

    bool operator==(const PeptideHit::PeakAnnotation& other) const;

    static void writePeakAnnotationsString_(String& annotation_string, std::vector<PeptideHit::PeakAnnotation> annotations);

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

    /// set information on (search engine) sub scores associated with this PSM (only used by pepXML)
    void setAnalysisResults(const std::vector<PepXMLAnalysisResult>& aresult);

    /// add information on (search engine) sub scores associated with this PSM (only used by pepXML)
    void addAnalysisResults(const PepXMLAnalysisResult& aresult);

    /// returns information on (search engine) sub scores associated with this PSM (only used by pepXML)
    std::vector<PepXMLAnalysisResult> getAnalysisResults() const;

    /// returns the PSM rank
    UInt getRank() const;

    /// sets the PSM rank (0 = top hit)
    void setRank(UInt newrank);

    /// returns the fragment annotations
    std::vector<PeptideHit::PeakAnnotation>& getPeakAnnotations();
    const std::vector<PeptideHit::PeakAnnotation>& getPeakAnnotations() const;


    /// sets the fragment annotations
    void setPeakAnnotations(std::vector<PeptideHit::PeakAnnotation> frag_annotations);

    /**
     * @brief Returns true if this hit is annotated as mapping to decoy sequences only. 
     *
     * Note: an unknown/unannotated state will yield false.
     */
    bool isDecoy() const;

    /** @brief Sets the target/decoy type for this peptide hit
     *
     * This method provides a type-safe way to annotate peptide hits with their
     * target/decoy status. Use TARGET_DECOY for peptides that match both target
     * and decoy protein sequences (these are treated as targets in FDR calculations).
     * Note: UNKNOWN should only be used in special cases where the status needs to
     * be explicitly marked as unknown.
     *
     * @param type The target/decoy classification:
     *   - TARGET: Only matches target proteins
     *   - DECOY: Only matches decoy proteins
     *   - TARGET_DECOY: Matches both target and decoy proteins
     *   - UNKNOWN: Target/decoy status is unknown (explicit unknown state)
     */
    void setTargetDecoyType(TargetDecoyType type);

    /** @brief Returns the target/decoy type for this peptide hit
     *
     * This method performs case-insensitive parsing of the "target_decoy" meta value
     * and returns the corresponding enum value. Returns UNKNOWN if the meta value
     * does not exist.
     *
     * @return The target/decoy classification:
     *   - TARGET: Only matches target proteins
     *   - DECOY: Only matches decoy proteins
     *   - TARGET_DECOY: Matches both target and decoy proteins
     *   - UNKNOWN: Target/decoy status not set (meta value missing)
     */
    TargetDecoyType getTargetDecoyType() const;

    //@}

    /// extracts the set of non-empty protein accessions from peptide evidences
    std::set<String> extractProteinAccessionsSet() const;

protected:
    AASequence sequence_;

    /// the score of the peptide hit
    double score_{};

    /// the charge of the peptide
    Int charge_{};

    /// information on the potential peptides observed through this PSM.
    std::vector<PeptideEvidence> peptide_evidences_;

    /// annotations of fragments in the corresponding spectrum
    std::vector<PeptideHit::PeakAnnotation> fragment_annotations_;

private:
    /// Get the number of analysis results stored as meta values (only for pepXML results)
    size_t getNumberOfAnalysisResultsFromMetaValues_() const;

    /// Extract analysis results from meta values (only for pepXML results)
    std::vector<PepXMLAnalysisResult> extractAnalysisResultsFromMetaValues_() const;
  };

  /// Stream operator
  OPENMS_DLLAPI std::ostream& operator<< (std::ostream& stream, const PeptideHit& hit);
} // namespace OpenMS

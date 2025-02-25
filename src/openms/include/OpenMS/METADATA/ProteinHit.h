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
#include <functional>
#include <set>
#include <map>

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
    @brief Representation of a protein hit

    It contains the fields score, score_type, rank, accession,
    sequence and coverage.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI ProteinHit :
    public MetaInfoInterface
  {
public:
    static const double COVERAGE_UNKNOWN; // == -1

    /// @name Hashes for ProteinHit
    //@{
    /// Hash of a ProteinHit based on its accession only!
    class OPENMS_DLLAPI ProteinHitAccessionHash
    {
    public:
      size_t operator()(const ProteinHit & p)
      {
        return std::hash<std::string>{}(p.getAccession());
      }

    };
    class OPENMS_DLLAPI ProteinHitPtrAccessionHash
    {
    public:
      size_t operator()(const ProteinHit * p)
      {
        return std::hash<std::string>{}(p->getAccession());
      }

    };
    //@}

    /// @name Comparators ProteinHit
    //@{
    /// Greater predicate for scores of hits
    class OPENMS_DLLAPI ScoreMore
    {
public:
      template <typename Arg>
      bool operator()(const Arg & a, const Arg & b)
      {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
        if (a.getScore() != b.getScore())
#pragma clang diagnostic pop
        {
          return a.getScore() > b.getScore();
        }
        return a.getAccession() > b.getAccession();
      }

    };

    /// Lesser predicate for scores of hits
    class OPENMS_DLLAPI ScoreLess
    {
public:
      template <typename Arg>
      bool operator()(const Arg & a, const Arg & b)
      {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
        if (a.getScore() != b.getScore())
#pragma clang diagnostic pop
        {
          return a.getScore() < b.getScore();
        }
        return a.getAccession() < b.getAccession();
      }

    };
    //@}

    /** @name Constructors and Destructor */
    //@{

    /// Default constructor
    ProteinHit();

    /// Values constructor
    ProteinHit(double score, UInt rank, String accession, String sequence);

    /// Copy constructor
    ProteinHit(const ProteinHit &) = default;

    /// Move constructor
    ProteinHit(ProteinHit&&) = default;

    //@}

    /// Assignment operator
    ProteinHit & operator=(const ProteinHit &) = default;

    /// Move assignment operator
    ProteinHit& operator=(ProteinHit&&) = default; // TODO: add noexcept (gcc 4.8 bug)

    /// Assignment for MetaInfo
    ProteinHit & operator=(const MetaInfoInterface & source);

    /// Equality operator
    bool operator==(const ProteinHit & rhs) const;

    /// Inequality operator
    bool operator!=(const ProteinHit & rhs) const;


    /** @name Accessors */
    //@{

    /// returns the score of the protein hit
    double getScore() const;

    /// returns the rank of the protein hit
    UInt getRank() const;

    /// returns the protein sequence
    const String & getSequence() const;

    /// returns the accession of the protein
    const String & getAccession() const;
    
    /// returns the description of the protein
    String getDescription() const;

    /// returns the coverage (in percent) of the protein hit based upon matched peptides
    double getCoverage() const;

    /// sets the score of the protein hit
    void setScore(const double score);

    /// sets the rank
    void setRank(UInt newrank);

    /// sets the protein sequence
    void setSequence(const String & sequence);
    void setSequence(String && sequence);

    /// sets the accession of the protein
    void setAccession(const String & accession);

    /// sets the description of the protein
    void setDescription(const String & description);

    /// sets the coverage (in percent) of the protein hit based upon matched peptides
    void setCoverage(const double coverage);

    /// returns the set of modified protein positions
    const std::set<std::pair<Size, ResidueModification> >& getModifications() const;

    /// sets the set of modified protein positions
    void setModifications(std::set<std::pair<Size, ResidueModification> >& mods);
    //@}

protected:
    double score_;       ///< the score of the protein hit
    UInt rank_;          ///< the position(rank) where the hit appeared in the hit list
    String accession_;   ///< the protein identifier
    String sequence_;    ///< the amino acid sequence of the protein hit
    double coverage_;    ///< coverage of the protein based upon the matched peptide sequences
    std::set<std::pair<Size, ResidueModification> > modifications_; ///< modified positions in a protein
  };

  /// Stream operator
  OPENMS_DLLAPI std::ostream& operator<< (std::ostream& stream, const ProteinHit& hit);

} // namespace OpenMS


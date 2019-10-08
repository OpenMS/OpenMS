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
    double score_;        ///< the score of the protein hit
    UInt rank_;          ///< the position(rank) where the hit appeared in the hit list
    String accession_;   ///< the protein identifier
    String sequence_;    ///< the amino acid sequence of the protein hit
    double coverage_;    ///< coverage of the protein based upon the matched peptide sequences
    std::set<std::pair<Size, ResidueModification> > modifications_; ///< modified positions in a protein
  };

} // namespace OpenMS


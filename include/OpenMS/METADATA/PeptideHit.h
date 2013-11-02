// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
      bool operator()(const Arg & a, const Arg & b)
      {
        return a.getScore() > b.getScore();
      }

    };

    /// Lesser predicate for scores of hits
    class OPENMS_DLLAPI ScoreLess
    {
public:
      template <typename Arg>
      bool operator()(const Arg & a, const Arg & b)
      {
        return a.getScore() < b.getScore();
      }

    };
    //@}

    /** @name Constructors and Destructor */
    //@{
    /// default constructor
    PeptideHit();

    /// values constructor
    PeptideHit(DoubleReal score,
               UInt rank,
               Int charge,
               const AASequence & sequence);

    /// copy constructor
    PeptideHit(const PeptideHit & source);

    /// destructor
    virtual ~PeptideHit();
    //@}

    /// assignment operator
    PeptideHit & operator=(const PeptideHit & source);

    /// Equality operator
    bool operator==(const PeptideHit & rhs) const;

    /// Inequality operator
    bool operator!=(const PeptideHit & rhs) const;

    /**	@name Accessors
    */
    //@{
    /// returns the score of the peptide hit
    DoubleReal getScore() const;

    /// returns the rank of the peptide hit
    UInt getRank() const;

    /// returns the peptide sequence without trailing or following spaces
    const AASequence & getSequence() const;

    /// returns the charge of the peptide
    Int getCharge() const;

    /// returns the corresponding protein accessions
    const std::vector<String> & getProteinAccessions() const;

    /// sets the corresponding protein accessions
    void setProteinAccessions(const std::vector<String> & accessions);

    /// sets the score of the peptide hit
    void setScore(DoubleReal score);

    /// sets the rank
    void setRank(UInt newrank);

    /// sets the peptide sequence
    void setSequence(const AASequence & sequence);

    /// sets the charge of the peptide
    void setCharge(Int charge);

    /// adds an accession of a protein which contains this peptide hit
    void addProteinAccession(const String & accession);

    /// sets the amino acid before the sequence
    void setAABefore(char acid);
    /// returns the amino acid before the sequence
    char getAABefore() const;

    /// sets the amino acid after the sequence
    void setAAAfter(char acid);
    /// returns the amino acid after the sequence
    char getAAAfter() const;


    //@}


protected:
    DoubleReal score_;          ///< the score of the peptide hit
    UInt rank_;                         ///< the position(rank) where the hit appeared in the hit list
    Int charge_;         ///< the charge of the peptide
    AASequence sequence_;                               ///< the amino acid sequence of the peptide hit
    char aa_before_;     ///< Amino acid before the sequence
    char aa_after_;     ///< Amino acid after the sequence
    std::vector<String> corresponding_protein_accessions_;     ///< the accessions of the corresponding proteins

  };

} // namespace OpenMS

#endif // OPENMS_METADATA_PEPTIDEHIT_H

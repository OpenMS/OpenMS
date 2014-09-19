// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

    class PeptideEvidence
    {
      public:
        static const int UNKNOWN_POSITION = -1;
        // Note: we use 0 as position of the N-terminus while e.g. mzTab or other formats start counting at 1.
        static const int N_TERMINAL_POSITION = 0;
        static const char UNKNOWN_AA = " ";
        // Note: we use '-' as in mzTab specification
        static const char N_TERMINAL_AA = "-";
        static const char C_TERMINAL_AA = "-";

        PeptideEvidence() :
          accession_(""),
          first_(UNKNOWN_POSITION),
          last_(UNKNOWN_POSITION),
          aa_before_(UNKNOWN_AA),
          aa_after_(UNKNOWN_AA)
        {
        }

        /// set the protein accession the peptide matches to. If not available set to empty string.
        void setAccession(String a);

        /// get the protein accession the peptide matches to. If not available the empty string is returned.
        String getAccession();

        /// set the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set to UNKNOWN_POSITION. N-terminal positions must be marked with N_TERMINAL_AA
        void setFirst(Int a);

        /// get the position in the protein (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned.
        Int getFirst() const;

        /// set the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set UNKNOWN_POSITION. C-terminal positions must be marked with C_TERMINAL_AA
        void setLast(Int a);

        /// get the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned.
        Int getLast() const;

        /// sets the amino acid single letter code before the sequence (preceeding amino acid in the protein). If not available, set to UNKNOWN_AA. If N-terminal set to N_TERMINAL_AA
        void setAABefore(char acid);

        /// returns the amino acid single letter code before the sequence (preceeding amino acid in the protein). If not available, UNKNOWN_AA is returned. If N-terminal, N_TERMINAL_AA is returned.
        char getAABefore() const;

        /// sets the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, set to UNKNOWN_AA. If C-terminal set to C_TERMINAL_AA
        void setAAAfter(char acid);

        /// returns the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, UNKNOWN_AA is returned. If C-terminal, C_TERMINAL_AA is returned.
        char getAAAfter() const;

     private:
        String accession_;
        Int first_;
        Int last_;
        char aa_before_;
        char aa_after_;
    };

    /** @name Constructors and Destructor */
    //@{
    /// default constructor
    PeptideHit();

    /// values constructor
    PeptideHit(double score,
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
    /// returns the peptide sequence without trailing or following spaces
    const AASequence & getSequence() const;

    /// sets the peptide sequence
    void setSequence(const AASequence & sequence);

    /// returns the charge of the peptide
    Int getCharge() const;

    /// sets the charge of the peptide
    void setCharge(Int charge);

    /// returns information on peptides (potentially) identified by this PSM
    const std::vector<String> & getPeptideEvidence() const;

    /// set information on peptides (potentially) identified by this PSM
    void setPeptideEvidence(const std::vector<String> & accessions);

    /// adds information on a peptide that is (potentially) identified by this PSM
    void addPeptideEvidence(const PeptideEvidence & peptide_evidence);

    /// returns the PSM score
    double getScore() const;

    /// sets the PSM score
    void setScore(double score);

    /// returns the PSM rank
    UInt getRank() const;

    /// sets the PSM rank
    void setRank(UInt newrank);
    //@}

protected:
    /// the score of the peptide hit
    double score_;

    /// the position(rank) where the hit appeared in the hit list
    UInt rank_;

    /// the charge of the peptide
    Int charge_;

    /// the amino acid sequence of the peptide hit
    AASequence sequence_;

    /// information on the potential peptides observed through this PSM. The name is borrowed from the mzIdent specification.
    std::vector<PeptideEvidence> peptide_evidences_;
  };

} // namespace OpenMS

#endif // OPENMS_METADATA_PEPTIDEHIT_H

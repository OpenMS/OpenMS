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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PEPTIDEEVIDENCE_H
#define OPENMS_METADATA_PEPTIDEEVIDENCE_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

/**
  @brief Representation of a peptide evidence.

  A peptide evidence object describes a single peptide to protein match.

  @ingroup Metadata
*/
  class OPENMS_DLLAPI PeptideEvidence
  {
public:
    static const int UNKNOWN_POSITION; // == -1

    // Note: we use 0 as position of the N-terminus while e.g. mzTab or other formats start counting at 1.
    static const int N_TERMINAL_POSITION;
    static const char UNKNOWN_AA; // PeptideEvidence::UNKNOWN_AA = 'X';

    // Note: we use '[' and ']' instead of  e.g. '-' as in mzTab specification
    static const char N_TERMINAL_AA;
    static const char C_TERMINAL_AA;

    /// constructor
    PeptideEvidence();

    /// constructor
    PeptideEvidence(const String& accession, Int start, Int end, char aa_before, char aa_after);

    /// copy constructor
    PeptideEvidence(const PeptideEvidence& source);

    /// destructor
    ~PeptideEvidence() {}
    //@}

    /// assignment operator
    PeptideEvidence& operator=(const PeptideEvidence& source);

    /// Equality operator
    bool operator==(const PeptideEvidence& rhs) const;
    
    /// Less operator
    bool operator<(const PeptideEvidence& rhs) const;

    /// not equal
    bool operator!=(const PeptideEvidence& rhs) const;

    /// start and end numbers in evidence represent actual numeric indices
    bool hasValidLimits() const;

    /// get the protein accession the peptide matches to. If not available the empty string is returned.
    const String& getProteinAccession() const;

    /// set the protein accession the peptide matches to. If not available set to empty string.
    void setProteinAccession(const String& s);

    /// set the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set to UNKNOWN_POSITION. N-terminal positions must be marked with N_TERMINAL_AA
    void setStart(const Int a);

    /// get the position in the protein (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned.
    Int getStart() const;

    /// set the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set UNKNOWN_POSITION. C-terminal positions must be marked with C_TERMINAL_AA
    void setEnd(const Int a);

    /// get the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned.
    Int getEnd() const;

    /// sets the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, set to UNKNOWN_AA. If N-terminal set to N_TERMINAL_AA
    void setAABefore(const char acid);

    /// returns the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, UNKNOWN_AA is returned. If N-terminal, N_TERMINAL_AA is returned.
    char getAABefore() const;

    /// sets the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, set to UNKNOWN_AA. If C-terminal set to C_TERMINAL_AA
    void setAAAfter(const char acid);

    /// returns the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, UNKNOWN_AA is returned. If C-terminal, C_TERMINAL_AA is returned.
    char getAAAfter() const;

protected:
    String accession_;

    Int start_;

    Int end_;

    char aa_before_;

    char aa_after_;
  };

}

#endif // OPENMS_METADATA_PEPTIDEEVIDENCE_H

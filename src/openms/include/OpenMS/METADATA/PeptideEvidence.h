// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

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

    /// Constructor
    PeptideEvidence();

    /// Constructor
    explicit PeptideEvidence(const String& accession, Int start=UNKNOWN_POSITION, Int end=UNKNOWN_POSITION, char aa_before=UNKNOWN_AA, char aa_after=UNKNOWN_AA);

    /// Copy constructor
    PeptideEvidence(const PeptideEvidence&) = default;

    /// Move constructor
    PeptideEvidence(PeptideEvidence&&) noexcept = default;

    /// Destructor
    ~PeptideEvidence() {}
    //@}

    /// Assignment operator
    PeptideEvidence& operator=(const PeptideEvidence&) = default;
    /// Move assignment operator
    PeptideEvidence& operator=(PeptideEvidence&&) = default; // TODO: add noexcept (gcc 4.8 bug)

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


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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PEPTIDEEVIDENCE_H
#define OPENMS_METADATA_PEPTIDEEVIDENCE_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  class AASequence;

  /**
    @brief Representation of a MzIdentML PeptideEvidence

        It contains information about the protein the peptide comes from
        and additional information were it is located in the protein.

        @todo implement CVParam

        @ingroup Metadata
  */
  class OPENMS_DLLAPI PeptideEvidence :
    public MetaInfoInterface
  {
public:

    /** @name Constructors and Destructor */
    //@{
    /// default constructor
    PeptideEvidence();

    /// values constructor
    PeptideEvidence(DoubleReal score,
                    UInt rank,
                    Int charge,
                    const AASequence & sequence);

    /// copy constructor
    PeptideEvidence(const PeptideEvidence & source);

    /// destructor
    virtual ~PeptideEvidence();
    //@}

    /// assignment operator
    PeptideEvidence & operator=(const PeptideEvidence & source);

    /// Equality operator
    bool operator==(const PeptideEvidence & rhs) const;

    /// Inequality operator
    bool operator!=(const PeptideEvidence & rhs) const;

    /**	@name Accessors
    */
    //@{
    /// returns the corresponding protein db sequence ref
    const String & getDBSequenceRef() const;

    /// sets the corresponding protein db sequence ref
    void setDBSequenceRef(const String & rhs);

    /// returns the translation table reference
    const String & getTranslationTableRef() const;

    /// sets the translation table reference
    void setTranslationTableRef(const String & rhs);

    /// start position in the sequence (xsd:int)
    void setStart(Int start);

    /// returns the start position in the sequence (xsd:int)
    Int getStart() const;

    /// sets the end position in the sequence
    void setEnd(Int end);

    /// returns the end position in the sequence
    Int getEnd() const;

    /// sets the amino acid before the sequence, "-" if N-terminal, "?" if not applicable (e.g. de novo)
    void setPre(char rhs);

    /// returns the amino acid before the sequence
    char getPre() const;

    /// sets the amino acid after the sequence, "-" if C-terminal, "?" if not applicable (e.g. de novo)
    void setPost(char rhs);

    /// returns the amino acid after the sequence
    char getPost() const;

    /// unique id of the file, set of files or repository
    void setId(const String & id);

    /// returns the unqiue id of the instance
    const String & getId() const;

    /// sets the potentially ambigous but human readable name
    void setName(const String & name);

    /// returns the human readable name
    const String & getName() const;

    /// sets the number of missed cleavages
    void setMissedCleavages(Int rhs);

    /// returns the number of missed cleavages
    Int getMissedCleavages() const;

    /// sets whether the hit is a decoy hit
    void setIsDecoy(bool is_decoy);

    /// returns the whether the hit is decoy
    bool getIsDecoy() const;

    /// Frame of the DB, e.g. from nucleic acids
    void setFrame(Int frame);

    /// returns the frame of the peptide evidence
    Int getFrame() const;
    //@}


protected:

    String db_sequence_ref_;
    String translation_table_ref_;
    Int start_;
    Int end_;
    char pre_;
    char post_;
    String id_;
    String name_;
    Int missed_cleavages_;
    bool is_decoy_;
    Int frame_;

  };

} // namespace OpenMS

#endif // OPENMS_METADATA_PEPTIDEEVIDENCE_H

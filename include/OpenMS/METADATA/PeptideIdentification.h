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

#ifndef OPENMS_METADATA_PEPTIDEIDENTIFICATION_H
#define OPENMS_METADATA_PEPTIDEIDENTIFICATION_H

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/ProteinHit.h>

#include <string>
#include <map>

namespace OpenMS
{
  /**
    @brief Represents the peptide hits for a spectrum

      This class is closely related to ProteinIdentification, which stores the protein hits
      and the general information about the identification run. More than one PeptideIdentification
      can belong to one ProteinIdentification. The general information about a
      PeptideIdentification has to be looked up in the correpsonding ProteinIndentification, using
      the unique <i>identifier</i> that links the two.

      When loading PeptideHit instances from a File, the retention time and mass-to-charge ratio
      of the precursor spectrum is stored in the MetaInfoInterface using the names 'MZ' and 'RT'.
      This information can be used to map the peptide hits to an MSExperiment, a FeatureMap
      or a ConsensusMap using the IDMapper class.

        @ingroup Metadata
  */
  class OPENMS_DLLAPI PeptideIdentification :
    public MetaInfoInterface
  {
public:

    ///Hit type definition
    typedef PeptideHit HitType;

    /// @name constructors,destructors,assignment operator
    //@{
    /// default constructor
    PeptideIdentification();
    /// destructor
    virtual ~PeptideIdentification();
    /// copy constructor
    PeptideIdentification(const PeptideIdentification & source);
    /// assignment operator
    PeptideIdentification & operator=(const PeptideIdentification & source);
    /// Equality operator
    bool operator==(const PeptideIdentification & rhs) const;
    /// Inequality operator
    bool operator!=(const PeptideIdentification & rhs) const;
    //@}

    /// returns the peptide hits
    const std::vector<PeptideHit> & getHits() const;
    /// Appends a peptide hit
    void insertHit(const PeptideHit & hit);
    /// Sets the peptide hits
    void setHits(const std::vector<PeptideHit> & hits);

    /// returns the peptide significance threshold value
    DoubleReal getSignificanceThreshold() const;
    /// setting of the peptide significance threshold value
    void setSignificanceThreshold(DoubleReal value);

    /// returns the peptide score type
    String getScoreType() const;
    /// sets the peptide score type
    void setScoreType(const String & type);

    /// returns the peptide score orientation
    bool isHigherScoreBetter() const;
    /// sets the peptide score orientation
    void setHigherScoreBetter(bool value);

    /// returns the identifier
    const String & getIdentifier() const;
    /// sets the indentifier
    void setIdentifier(const String & id);

    /// Sorts the hits by score and assigns ranks coording to the scores
    void assignRanks();

    /**
         @brief Sorts the hits by score

         Sorting takes the score orientation (@p higher_score_better_) into account, i.e. after sorting, the best-scoring hit is the first.
    */
    void sort();

    /// Returns if this PeptideIdentification result is empty
    bool empty() const;

    ///@name Methods for linking peptide and protein hits
    //@{

    /// finds and inserts all peptide hits which reference to a given protein accession
    void getReferencingHits(const String & protein_accession, std::vector<PeptideHit> & peptide_hits) const;

    /// finds and inserts all peptide hits which reference to a given list of protein accessions
    void getReferencingHits(const std::vector<String> & accessions, std::vector<PeptideHit> & peptide_hits) const;

    /// finds and inserts all peptide hits which reference to a given list of proteins (via their accessions)
    void getReferencingHits(const std::vector<ProteinHit> & protein_hits, std::vector<PeptideHit> & peptide_hits) const;

    /// the complement of the above
    void getNonReferencingHits(const String & protein_accession, std::vector<PeptideHit> & peptide_hits) const;
    void getNonReferencingHits(const std::vector<String> & accessions, std::vector<PeptideHit> & peptide_hits) const;
    void getNonReferencingHits(const std::vector<ProteinHit> & protein_hits, std::vector<PeptideHit> & peptide_hits) const;
    //@}

protected:
    String id_;                                                          ///< Identifier by which ProteinIdentification and PeptideIdentification are matched
    std::vector<PeptideHit> hits_;               ///< A list containing the peptide hits
    DoubleReal significance_threshold_;              ///< the peptide significance threshold
    String score_type_;                                      ///< The score type (Mascot, Sequest, e-value, p-value)
    bool higher_score_better_;                       ///< The score orientation
  };

} //namespace OpenMS
#endif // OPENMS_METADATA_PEPTIDEIDENTIFICATION_H

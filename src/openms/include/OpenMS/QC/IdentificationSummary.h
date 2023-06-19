// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <vector>
/**
 * @brief Detected Proteins/Peptides as a Proteomics QC metric
 *
 * Simple class to return a summary of detected proteins/peptides
 * from a given idXML file.
 *
 */

namespace OpenMS
{
  class OPENMS_DLLAPI IdentificationSummary : public QCBase
  {
  public:
    /// Constructor
    IdentificationSummary() = default;

    /// Destructor
    virtual ~IdentificationSummary() = default;

    // small struct for unique peptide / protein identifications (considering sequence only)
    // count: number of unique identifications, fdr_threshold: significance threshold if score type is FDR, else -1
    struct OPENMS_DLLAPI UniqueID {
      UInt count = 0;
      float fdr_threshold = -1.0;
    };

    // stores identification summary values calculated by compute function
    struct OPENMS_DLLAPI Result {
      UInt peptide_spectrum_matches = 0;
      UniqueID unique_peptides;
      UniqueID unique_proteins;
      float missed_cleavages_mean = 0;
      double protein_hit_scores_mean = 0;
      double peptide_length_mean = 0;

      bool operator==(const Result& rhs) const;
    };

    /**
   @brief computes a summary of an idXML file

   @param prot_ids vector with ProteinIdentifications
   @param pep_ids vector with PeptideIdentifications
   @return result object with summary values:
           total number of PSM (peptide_spectrum_matches),
           number of identified peptides with given FDR threshold (unique_peptides),
           number of identified proteins with given FDR threshold (unique_proteins),
           missed cleavages mean (missed_cleavages_mean),
           identification score mean of protein hits (protein_hit_scores_mean),
           identified peptide lengths mean (peptide_length_mean)

   **/
    Result compute(std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);

    const String& getName() const override;

    QCBase::Status requirements() const override;

  private:
    const String name_ = "Summary of detected Proteins and Peptides from idXML file";
  };
} // namespace OpenMS

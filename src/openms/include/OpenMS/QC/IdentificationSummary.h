// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

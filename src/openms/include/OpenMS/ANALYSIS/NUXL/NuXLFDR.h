// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLReport.h>
#include <vector>

namespace OpenMS
{

/// @brief adapted FDR calculation for NA cross-links
class OPENMS_DLLAPI NuXLFDR
{
  public:
    explicit NuXLFDR(size_t report_top_hits);

    // split by meta value "NuXL:isXL" == 0
    void splitIntoPeptidesAndXLs(const PeptideIdentificationList& peptide_ids, 
      PeptideIdentificationList& pep_pi, 
      PeptideIdentificationList& xl_pi) const;

    void mergePeptidesAndXLs(const PeptideIdentificationList& pep_pi, 
      const PeptideIdentificationList& xl_pi, 
      PeptideIdentificationList& peptide_ids) const;

    // calculate PSM-level q-values (irrespective of XL/non-XL class)
    void QValueAtPSMLevel(PeptideIdentificationList& peptide_ids) const;

    // calculate PSM-level q-values for XL and non-XL class separately.
    void calculatePeptideAndXLQValueAtPSMLevel(const PeptideIdentificationList& peptide_ids, 
      PeptideIdentificationList& pep_pi, 
      PeptideIdentificationList& xl_pi) const;

    // calculate separate FDRs, filter decoys, write PSM and protein reports
    void calculatePeptideAndXLQValueAndFilterAtPSMLevel(
      const std::vector<ProteinIdentification>& protein_ids,
      const PeptideIdentificationList& peptide_ids, 
      PeptideIdentificationList& pep,
      double peptide_PSM_qvalue_threshold,
      double peptide_peptide_qvalue_threshold,
      PeptideIdentificationList& xl_pi,
      std::vector<double> xl_PSM_qvalue_thresholds,
      std::vector<double> xl_peptidelevel_qvalue_thresholds,
      const String& out_idxml,
      int decoy_factor) const;

  private:
    size_t report_top_hits_;
};

}



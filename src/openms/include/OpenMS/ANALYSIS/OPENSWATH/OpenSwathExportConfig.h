// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportData.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/config.h>

namespace OpenMS
{
  /**
    @brief Shared filter configuration for OpenSWATH export readers.

    Defaults mirror the corresponding PyProphet export CLI defaults

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathExportFilterConfig
  {
    OpenSwathIPFExportMode ipf_mode = OpenSwathIPFExportMode::Peptidoform;
    /// Report aggregated transition-level quantification
    bool transition_quantification = true;
    /// Maximum PEP to retain scored transitions for quantification (requires transition-level scoring)
    double max_transition_pep = 0.7;
    /// Filter results to maximum run-specific peptidoform-level PEP
    double ipf_max_peptidoform_pep = 0.4;
    /// Filter results to maximum run-specific peak group-level q-value
    double max_rs_peakgroup_qvalue = 0.05;
    
    bool peptide = true;
    double max_global_peptide_qvalue = 0.01;
    bool protein = true;
    double max_global_protein_qvalue = 0.01;
    bool gene = false;
    double max_global_gene_qvalue = 0.01;
    bool use_alignment = false;
    double max_alignment_pep = 0.7;
    bool exclude_decoys = true;
  };

  /**
    @brief Configuration for user-facing feature results exports.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathResultsExportConfig
  {
    OpenSwathExportFilterConfig filters;
    OpenSwathExportFileFormat format = OpenSwathExportFileFormat::TSV;
  };

  /**
    @brief Configuration for matrix exports.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathMatrixExportConfig
  {
    OpenSwathExportFilterConfig filters;
    OpenSwathMatrixLevel level = OpenSwathMatrixLevel::Peptide;
    OpenSwathMatrixNormalization normalization = OpenSwathMatrixNormalization::None;
    Size top_n = 3;
    bool consistent_top = true;
    OpenSwathExportFileFormat format = OpenSwathExportFileFormat::TSV;
  };

  /**
    @brief Configuration for scored feature Parquet exports.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathParquetExportConfig
  {
    OpenSwathExportFilterConfig filters;
    bool include_transition_data = true;
  };
} // namespace OpenMS

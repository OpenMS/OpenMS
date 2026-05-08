// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathInferenceData.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/config.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Configuration for OpenSwath analyte level inference global and local error estimation.

    Defaults mirror the current PyProphet CLI defaults.

    @ingroup OpenSwath
  */
  struct OPENMS_DLLAPI ErrorEstimationConfig
  {
    bool parametric = false;
    bool pfdr = false;
    std::vector<double> pi0_lambda {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45};
    String pi0_method = "bootstrap";
    Int32 pi0_smooth_df = 3;
    bool pi0_smooth_log_pi0 = false;
    bool lfdr_truncate = true;
    bool lfdr_monotone = true;
    String lfdr_transformation = "probit";
    double lfdr_adj = 1.5;
    double lfdr_eps = 1e-8;
  };

  /**
    @brief Configuration for peptidoform inference.

    Defaults mirror the current PyProphet CLI defaults.

    @ingroup OpenSwath
  */
  struct OPENMS_DLLAPI PeptidoformInferenceConfig
  {
    bool ipf_ms1_scoring = true;
    bool ipf_ms2_scoring = true;
    bool ipf_h0 = true;
    bool ipf_grouped_fdr = false;
    double ipf_max_peakgroup_pep = 0.7;
    double ipf_max_precursor_pep = 0.7;
    double ipf_max_precursor_peakgroup_pep = 0.4;
    double ipf_max_transition_pep = 0.6;
    bool propagate_signal_across_runs = false;
    double ipf_min_alignment_mapping_confidence = 0.5;
    double across_run_confidence_threshold = 0.5;
  };

  /**
    @brief Configuration for peptide-, protein-, and gene-level context inference.

    @ingroup OpenSwath
  */
  struct OPENMS_DLLAPI LevelContextInferenceConfig
  {
    InferenceLevel level = InferenceLevel::Peptide;
    InferenceContext context = InferenceContext::RunSpecific;
    ErrorEstimationConfig error;
  };
} // namespace OpenMS

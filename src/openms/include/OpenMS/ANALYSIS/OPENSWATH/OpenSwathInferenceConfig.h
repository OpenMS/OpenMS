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

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI ErrorEstimationConfig
  {
    /// If true, perform parametric estimation of p-values.
    bool parametric = false;
    /// If true, compute positive false discovery rate (pFDR) instead of FDR.
    bool pfdr = false;
    /// Lambda grid for pi0 estimation for non-parametric p-value estimation.
    std::vector<double> pi0_lambda {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45};
    /// Method for choosing the pi0 tuning parameter ("smoother", "bootstrap").
    String pi0_method = "bootstrap";
    /// Degrees of freedom for smoother-based pi0 estimation.
    Int32 pi0_smooth_df = 3;
    /// If true, smooth log(pi0) instead of pi0 directly.
    bool pi0_smooth_log_pi0 = false;
    /// If true, local FDR values above one are truncated to one.
    bool lfdr_truncate = true;
    /// If true, enforce non-decreasing local FDR values with increasing p-values.
    bool lfdr_monotone = true;
    /// Transformation applied to p-values for local FDR estimation ("probit", "logit").
    String lfdr_transformation = "probit";
    /// Smoothing bandwidth multiplier used during local FDR density estimation.
    double lfdr_adj = 1.5;
    /// Tail clipping threshold used during local FDR calculation.
    double lfdr_eps = 1e-8;
  };

  /**
    @brief Configuration for peptidoform inference.

    Defaults mirror the current PyProphet CLI defaults.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI PeptidoformInferenceConfig
  {
    /// If true, use MS1 precursor data for IPF.
    bool ipf_ms1_scoring = false;
    /// If true, use MS2 precursor data for IPF.
    bool ipf_ms2_scoring = false;
    /// If true, include the null hypothesis that no peptidoform explains the peak group.
    bool ipf_h0 = true;
    /// If true, compute grouped FDR instead of pooled FDR.
    bool ipf_grouped_fdr = false;
    /// Maximum peakgroup PEP to consider scored peak groups in IPF.
    double ipf_max_peakgroup_pep = 0.7;
    /// Maximum precursor PEP to consider scored precursor evidence in IPF.
    double ipf_max_precursor_pep = 0.7;
    /// Maximum precursor-peakgroup PEP retained for transition-level IPF.
    double ipf_max_precursor_peakgroup_pep = 0.4;
    /// Maximum transition PEP to consider scored transition evidence in IPF.
    double ipf_max_transition_pep = 0.6;
    /// If true, propagate confident transition evidence across aligned runs.
    bool propagate_signal_across_runs = false;
    /// Minimum mapping confidence retained from FEATURE_MS2_ALIGNMENT_CANDIDATE.
    double ipf_min_alignment_mapping_confidence = 0.5;
    /// Maximum transition PEP eligible for across-run evidence propagation.
    double across_run_confidence_threshold = 0.5;
  };

  /**
    @brief Configuration for peptide-, protein-, and gene-level context inference.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI LevelContextInferenceConfig
  {
    /// The analyte level to perform inference on.
    InferenceLevel level = InferenceLevel::Peptide;
    /// The context in which to perform inference (global, experiment-wide, or run-specific).
    InferenceContext context = InferenceContext::RunSpecific;
    /// Error-estimation settings applied during inference.
    ErrorEstimationConfig error;
  };
} // namespace OpenMS

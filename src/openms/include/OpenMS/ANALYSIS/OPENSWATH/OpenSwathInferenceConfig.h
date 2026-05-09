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

    @var bool parametric If true, perform parametric estimation of p-values.
    @var bool pfdr If true, compute positive false discovery rate (pFDR) instead of false discovery rate (FDR).
    @var std::vector<double> pi0_lambda Lambda grid for pi0 estimation for non-parametric estimation of p-values.
    @var String pi0_method Method for choosing tuning parameter in pi0 estimation ("smoother", "bootstrap").
    @var Int32 pi0_smooth_df Degrees of freedom for smoother-based pi0 estimation.
    @var bool pi0_smooth_log_pi0 If true, smooth log(pi0) instead of pi0 directly.
    @var bool lfdr_truncate If true, local FDR values > 1 are set to 1.
    @var bool lfdr_monotone If true, ensure local FDR values are non-decreasing with p-values.
    @var String lfdr_transformation Transformation applied to p-values for local FDR estimation ("probit", "logit").
    @var double lfdr_adj Smoothing bandwidth multiplier used in density estimation.
    @var double lfdr_eps Threshold for p-value tails in local FDR calculation.

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

    @var bool ipf_ms1_scoring If true, use MS1 precursor data for IPF.
    @var bool ipf_ms2_scoring If true, use MS2 precursor data for IPF.
    @var bool ipf_h0 If true, include possibility that peak groups are not covered by peptidoform space.
    @var bool ipf_grouped_fdr If true, compute grouped FDR instead of pooled FDR to better support data where peak groups are evaluated to originate from very heterogeneous numbers of peptidoforms.
    @var double ipf_max_peakgroup_pep Maximum PEP to consider scored peak groups in IPF.
    @var double ipf_max_precursor_pep Maximum PEP to consider scored precursors in IPF.
    @var double ipf_max_precursor_peakgroup_pep Maximum BHM layer 1 integrated precursor peakgroup PEP to consider in IPF.
    @var double ipf_max_transition_pep Maximum PEP to consider scored transitions in IPF.
    @var bool propagate_signal_across_runs If true, propagate confident transition evidence across aligned runs (requires aligned features cross runs).
    @var double ipf_min_alignment_mapping_confidence Minimum mapping confidence retained for across-run propagation from FEATURE_MS2_ALIGNMENT_CANDIDATE. Only alignment mappings with a mapping confidence above this threshold will be considered for across-run propagation of transition evidence. 
    @var double across_run_confidence_threshold Maximum PEP to consider for propagating signal across runs for aligned features.

    @ingroup OpenSwath
  */
  struct OPENMS_DLLAPI PeptidoformInferenceConfig
  {
    bool ipf_ms1_scoring = false;
    bool ipf_ms2_scoring = false;
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

      @var InferenceLevel level The analyte level to perform inference on.
      @var InferenceContext context The context in which to perform inference (global, experiment-wide, or run-specific).
      @var ErrorEstimationConfig error Configuration for error estimation to be performed as part of inference.

    @ingroup OpenSwath
  */
  struct OPENMS_DLLAPI LevelContextInferenceConfig
  {
    InferenceLevel level = InferenceLevel::Peptide;
    InferenceContext context = InferenceContext::RunSpecific;
    ErrorEstimationConfig error;
  };
} // namespace OpenMS

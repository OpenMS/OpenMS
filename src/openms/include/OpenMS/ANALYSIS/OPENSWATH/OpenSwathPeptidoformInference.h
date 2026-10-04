// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathInferenceConfig.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathInferenceData.h>
#include <OpenMS/config.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief OpenSwath inference of peptidoforms (IPF).

    This class is file-format agnostic and works only on compact typed rows.

    See: Rosenberger, G., Liu, Y., Röst, H. et al. Inference and quantification of peptidoforms in large sample cohorts by SWATH-MS. Nat Biotechnol 35, 781–788 (2017). https://doi.org/10.1038/nbt.3908

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathPeptidoformInference
  {
  public:
    /**
      @brief Compact Bayesian model row used by the public helper methods.

      The meaning of @p hypothesis depends on the inference layer:
      - precursor inference uses hypothesis @c 1 for the true precursor state and
        hypothesis @c 0 for the false precursor state
      - transition / peptidoform inference uses peptide IDs as hypotheses and
        reserves @c -1 for the H0/null hypothesis

      @ingroup TargetedQuantitation
    */
    struct BayesianModelRow
    {
      Int64 feature_id = -1;
      Int64 hypothesis = -1;
      double prior = 0.0;
      double evidence = 0.0;
    };

    /**
      @brief Posterior probability row returned by @ref applyBM.

      The @p hypothesis value follows the same conventions as in
      @ref BayesianModelRow.

      @ingroup TargetedQuantitation
    */
    struct PosteriorRow
    {
      Int64 feature_id = -1;
      Int64 hypothesis = -1;
      double posterior = 0.0;
    };

    /**
      @brief Compute model FDR/q-values from posterior error probabilities.

      This reuses OpenMS multiple-testing utilities and follows the PyProphet IPF
      convention of computing model FDR from sorted PEP values with maximum-rank
      tie handling.

      @param pep_values Posterior error probabilities for one model-FDR family
      @return Model FDR / q-values in the original input order
    */
    static std::vector<double> computeModelFDR(const std::vector<double>& pep_values);

    /**
      @brief Build precursor-layer Bayesian model rows.

      For each feature this emits rows for the true (@c 1) and false (@c 0)
      precursor hypotheses using the configured MS1/MS2 precursor evidence.

      @param rows Peakgroup and optional precursor-evidence rows
      @return Bayesian-model rows for precursor inference
    */
    static std::vector<BayesianModelRow> preparePrecursorBM(const std::vector<IPFPrecursorRow>& rows);

    /**
      @brief Apply a compact Bayesian model and return posterior probabilities.

      Rows are grouped by feature and hypothesis. Within each group, evidence
      values are multiplied and the smallest prior in the group is retained,
      matching the PyProphet implementation being ported.

      @param rows Bayesian-model rows
      @return Posterior probabilities per feature/hypothesis pair
    */
    static std::vector<PosteriorRow> applyBM(const std::vector<BayesianModelRow>& rows);

    /**
      @brief Perform precursor-layer inference and return precursor peakgroup PEPs.

      This step combines peakgroup evidence with optional MS1/MS2 precursor
      evidence. The returned precursor peakgroup PEP is then used as the prior
      for the transition/peptidoform layer.

      @param precursor_rows Peakgroup/precursor evidence rows
      @param config IPF configuration
      @return One precursor peakgroup PEP per retained feature
    */
    std::vector<IPFPrecursorProbabilityRow> precursorInference(const std::vector<IPFPrecursorRow>& precursor_rows,
                                                               const PeptidoformInferenceConfig& config) const;

    /**
      @brief Perform complete peptidoform inference.

      @param precursor_rows Peakgroup/precursor evidence
      @param transition_rows Transition-level evidence
      @param alignment_rows Optional alignment-group memberships for across-run propagation
      @param config IPF configuration

      @return Final IPF result rows without the H0 hypothesis
    */
    std::vector<IPFResultRow> infer(const std::vector<IPFPrecursorRow>& precursor_rows,
                                    const std::vector<IPFTransitionRow>& transition_rows,
                                    const std::vector<IPFAlignmentRow>& alignment_rows,
                                    const PeptidoformInferenceConfig& config) const;
  };
} // namespace OpenMS

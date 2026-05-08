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

    @ingroup OpenSwath
  */
  class OPENMS_DLLAPI PeptidoformInference
  {
  public:
    /**
      @brief Compact Bayesian model row used by the public helper methods.

      @ingroup OpenSwath
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

      @ingroup OpenSwath
    */
    struct PosteriorRow
    {
      Int64 feature_id = -1;
      Int64 hypothesis = -1;
      double posterior = 0.0;
    };

    /// Reuse OpenMS multiple-testing utilities for model FDR computation.
    static std::vector<double> computeModelFDR(const std::vector<double>& pep_values);

    /// Build precursor-layer Bayesian model rows.
    static std::vector<BayesianModelRow> preparePrecursorBM(const std::vector<IPFPrecursorRow>& rows);

    /// Apply the Bayesian model and return posterior probabilities per feature and hypothesis.
    static std::vector<PosteriorRow> applyBM(const std::vector<BayesianModelRow>& rows);

    /// Perform precursor-layer inference and return precursor peakgroup PEPs.
    std::vector<IPFPrecursorProbabilityRow> precursorInference(const std::vector<IPFPrecursorRow>& precursor_rows,
                                                               const PeptidoformInferenceConfig& config) const;

    /**
      @brief Perform complete peptidoform inference.

      @param precursor_rows Peakgroup/precursor evidence
      @param transition_rows Transition-level evidence
      @param config IPF configuration
      @return Final IPF result rows without the H0 hypothesis
    */
    std::vector<IPFResultRow> infer(const std::vector<IPFPrecursorRow>& precursor_rows,
                                    const std::vector<IPFTransitionRow>& transition_rows,
                                    const std::vector<IPFAlignmentRow>& alignment_rows,
                                    const PeptidoformInferenceConfig& config) const;
  };
} // namespace OpenMS

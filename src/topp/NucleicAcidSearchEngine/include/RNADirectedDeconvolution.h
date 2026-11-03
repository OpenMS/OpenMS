// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Sam Wein $
// $Authors: Sam Wein $
// --------------------------------------------------------------------------
#pragma once
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>

namespace OpenMS
{
  /// NASE-specific selection of precursor hypotheses and competing envelope evidence.
  class RNADirectedDeconvolution : public SpectralDeconvolution
  {
  public:
    /** @brief Deconvolve queries against complete parent MS1 evidence.
      Queries do not protect candidates from ordinary FLASH filtering.
      Empty queries run ordinary FLASH. Query state is cleared after every call.
      Competition uses the union of query/evidence intervals within the batch;
      subset-invariant callers must retain the same batch or use ordinary FLASH.
      @param[in] spectrum Complete parent MS1 spectrum.
      @param[in] scan_number Parent scan number.
      @param[in] queries Selected-ion m/z and signed charge; zero charge considers all charges.
      @throws Exception::IllegalArgument For nonempty non-MS1 queries or invalid m/z.
    */
    void deconvolve(const MSSpectrum& spectrum, int scan_number, const std::vector<std::pair<double, int>>& queries);

  private:
    std::vector<Precursor> queries_;
    std::optional<MassRanges> selectCandidateMassRanges_(const CandidateSelectionContext& context,
      const DeconvolvedSpectrum* evidence_groups, const std::vector<Size>& valid) const override;
  };
}

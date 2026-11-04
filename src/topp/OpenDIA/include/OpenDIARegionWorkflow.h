// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>

namespace OpenDIA
{
/** @brief OpenDIA extraction and scoring over shared, leased scan columns.
  Peak picking, scores and output conversion retain the shared workflow implementation.
  Evidence can be decoded again after eviction from the bounded raw cache.
  OSW candidate scoring is ordered by picked apex RT within the existing assay batch.
  Complete chromatograms are retained until that batch is scored, not across the run.
  Non-native accesses and non-tophat filters retain the shared extractor.
*/
class RegionWorkflow : public OpenMS::OpenSwathWorkflow
{
public:
  using OpenMS::OpenSwathWorkflow::OpenSwathWorkflow;

  /** @brief Extract one prepared acquisition lane using its exact mass snapshot.
    @param[in] input Lane access, optionally wrapped by one mass transformation.
    @param[in,out] output Existing chromatogram arrays, appended as in the shared extractor.
    @param[in] coordinates Prepared, m/z-sorted extraction coordinates.
    @param[in] parameters Physical extraction windows and filter.
  */
  void extractChromatograms_(const OpenSwath::SpectrumAccessPtr& input,
                             std::vector<OpenSwath::ChromatogramPtr>& output,
                             const std::vector<OpenMS::ChromatogramExtractor::ExtractionCoordinates>& coordinates,
                             const OpenMS::ChromExtractParams& parameters) override;

  /** @brief Reuse budgeted native extraction columns for read-only scoring.
    @param[in] input Final lane access, optionally wrapped by one mass transformation.
    @return Shared evidence access, or the unchanged input for unsupported readers.
  */
  OpenSwath::SpectrumAccessPtr spectrumAccessForScoring_(const OpenSwath::SpectrumAccessPtr& input) const override;

protected:
  /// Score picked candidates in physical RT order within each existing OSW batch.
  bool useCandidateRTScoring_() const override { return true; }
};
} // namespace OpenDIA

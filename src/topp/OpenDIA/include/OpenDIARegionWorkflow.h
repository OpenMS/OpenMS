// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>

namespace OpenDIA
{
/** @brief OpenDIA extraction over leased native columns.
  Scoring and output conversion retain the shared workflow implementation.
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
};
} // namespace OpenDIA

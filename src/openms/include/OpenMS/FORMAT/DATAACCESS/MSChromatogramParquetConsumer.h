// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <memory>

namespace OpenMS
{
  class MSChromatogramParquetConsumerImpl;

  /**
    @brief Writes chromatograms to a Parquet file with a PyProphet-compatible schema.

    The schema includes precursor/transition metadata, RT/intensity arrays and
    compression flags. Additional columns are run_id, source_file, and ms_level.
  */
  class OPENMS_DLLAPI MSChromatogramParquetConsumer : public Interfaces::IMSDataConsumer
  {
  public:
    MSChromatogramParquetConsumer(const String& filename,
                                  UInt64 run_id,
                                  const String& source_file,
                                  const OpenSwath::LightTargetedExperiment& transition_exp);

    ~MSChromatogramParquetConsumer() override;

    void consumeSpectrum(SpectrumType& s) override;
    void consumeChromatogram(ChromatogramType& c) override;
    void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) override;
    void setExperimentalSettings(const ExperimentalSettings& exp) override;

  private:
    std::unique_ptr<MSChromatogramParquetConsumerImpl> impl_;
  };
} // namespace OpenMS

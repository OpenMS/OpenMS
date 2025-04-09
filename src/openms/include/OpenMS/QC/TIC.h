// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/QC/QCBase.h>

/**
 * @brief Total Ion Count (TIC) as a QC metric
 *
 * Simple class to calculate the TIC of an MSExperiment.
 * Allows for multiple usage, because each calculated TIC is
 * stored internally. Those results can then be returned using
 * getResults().
 *
 */

namespace OpenMS
{
  class MzTabMetaData;
  class MSExperiment;
  class MSChromatogram;

  class OPENMS_DLLAPI TIC : public QCBase
  {
  public:
    /// Constructor
    TIC() = default;

    /// Destructor
    virtual ~TIC() = default;

    // stores TIC values calculated by compute function
    struct OPENMS_DLLAPI Result {
      std::vector<UInt> intensities; // TIC intensities
      std::vector<float> relative_intensities;
      std::vector<float> retention_times; // TIC RTs in seconds
      UInt area = 0;                      // Area under TIC
      UInt fall = 0;                      // MS1 signal fall (10x) count
      UInt jump = 0;                      // MS1 signal jump (10x) count

      bool operator==(const Result& rhs) const;
    };


    /**
    @brief Compute Total Ion Count and applies the resampling algorithm, if a bin size in RT seconds greater than 0 is given.

    All MS1 TICs within a bin are summed up.

    @param exp Peak map to compute the MS1 tick from
    @param bin_size RT bin size in seconds
    @param ms_level MS level of spectra for calculation
    @return result struct with with computed QC metrics: intensities, RTs (in seconds), area under TIC, 10x MS1 signal fall, 10x MS1 signal jump

    **/
    Result compute(const MSExperiment& exp, float bin_size = 0, UInt ms_level = 1);

    const String& getName() const override;

    const std::vector<MSChromatogram>& getResults() const;

    QCBase::Status requirements() const override;

    /// append QC data for given metrics to mzTab's MTD section
    void addMetaDataMetricsToMzTab(MzTabMetaData& meta, std::vector<Result>& tics);

  private:
    const String name_ = "TIC";
  };
} // namespace OpenMS

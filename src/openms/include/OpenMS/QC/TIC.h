// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/QC/QCBase.h>

namespace OpenMS
{
  class MzTabMetaData;
  class MSExperiment;
  class MSChromatogram;

  /**
    @brief Total Ion Current (TIC) as a QC metric.

    Builds a Total Ion Current trace at a chosen MS level from an
    @ref MSExperiment and returns a @ref Result that bundles the
    per-bin intensities, the relative intensities (0 -- 100), the
    retention times, the area under the TIC and two stability counters
    that flag abrupt 10-fold intensity changes between consecutive bins.
    Each call to @ref compute is independent; results are returned to
    the caller and not retained by the instance.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI TIC : public QCBase
  {
  public:
    /// Default constructor.
    TIC() = default;

    /// Destructor.
    virtual ~TIC() = default;

    /**
      @brief Computed Total Ion Current trace plus stability counters.

      Populated by @ref compute and otherwise default-constructed.
      Empty when the input @c MSExperiment had no spectra of the
      requested MS level.
    */
    struct OPENMS_DLLAPI Result {
      std::vector<UInt> intensities;       ///< Per-bin TIC intensities.
      std::vector<float> relative_intensities; ///< Per-bin intensities rescaled to [0, 100] of the maximum bin (0 when the maximum is 0).
      std::vector<float> retention_times;  ///< Per-bin retention times (seconds).
      UInt area = 0;                       ///< Sum of @ref intensities (area under the TIC).
      UInt fall = 0;                       ///< Number of consecutive bins where the intensity drops by at least 10x.
      UInt jump = 0;                       ///< Number of consecutive bins where the intensity rises by at least 10x.

      bool operator==(const Result& rhs) const;
    };


    /**
      @brief Compute the Total Ion Current trace and stability metrics.

      Builds the TIC at the requested MS level. When @p bin_size is
      positive, intensities are summed within @p bin_size second
      retention-time bins; when @p bin_size is @c 0, no resampling is
      applied and the trace contains one entry per source spectrum.

      The returned @ref Result is empty when @p exp has no spectra of
      MS level @p ms_level.

      @param[in] exp      Source experiment.
      @param[in] bin_size RT bin size in seconds; @c 0 disables binning.
      @param[in] ms_level MS level whose spectra are aggregated.
      @return Computed metrics. The @ref Result::intensities and
              @ref Result::retention_times are empty when no spectra
              of @p ms_level are present.
    */
    Result compute(const MSExperiment& exp, float bin_size = 0, UInt ms_level = 1);

    /// Name of this QC metric (@c "TIC").
    const String& getName() const override;

    const std::vector<MSChromatogram>& getResults() const;

    /// Input-data requirements of @ref compute (raw mzML data).
    QCBase::Status requirements() const override;

    /**
      @brief Append the TIC traces from @p tics to the MTD section of @p meta as custom parameters.

      For every non-empty entry in @p tics one custom parameter is added
      to @p meta, named @c TIC_1, @c TIC_2, ... (1-based index into
      @p tics), with the controlled-vocabulary label @c "total ion current"
      (@c MS:1000285). The value is a flat list of
      @c [rt0,int0,rt1,int1,...] pairs taken from @ref Result::retention_times
      and @ref Result::intensities.

      @param[in,out] meta mzTab metadata section to extend.
      @param[in]     tics TIC results to serialise; empty entries are skipped.
    */
    void addMetaDataMetricsToMzTab(MzTabMetaData& meta, const std::vector<Result>& tics);

  private:
    const String name_ = "TIC";
  };
} // namespace OpenMS

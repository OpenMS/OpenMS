// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/PeakMapExtractor.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS::OpenDIA
{
  /// One precursor whose raw MS1 and MS2 points are kept for inspection.
  struct PeakMapTarget
  {
    /// Compound id; precursor rows use the native id "<id>_Precursor_i<k>".
    std::string precursor_id;
    double precursor_mz = 0.0;
    /// Precursor charge; MS1 isotopes beyond the monoisotope need a positive charge.
    int charge = 0;
    /// Extraction centre in run seconds (the calibrated RT when a calibration is applied).
    double center_rt = std::numeric_limits<double>::quiet_NaN();
    bool has_rt = false;
    /// Precursor inverse reduced mobility (1/K0).
    double precursor_im = 0.0;
    bool has_im = false;
    /// Transition native id and extraction m/z, one MS2 peak map each.
    std::vector<std::pair<std::string, double>> fragments;
  };

  struct PeakMapOptions
  {
    /// Half width of every m/z extraction window, in ppm of the target m/z.
    double mz_tolerance_ppm = 20.0;
    /// Half width of the RT window around center_rt in seconds. Zero, or a target without RT, keeps the whole run.
    double rt_half_window_s = 0.0;
    /// Half width of the 1/K0 filter. Zero keeps every mobility (the mobility is still recorded).
    double im_half_window = 0.0;
    /// Additional MS1 isotopes after the monoisotope (k = 1..n, spaced 1.0033548 / charge).
    std::size_t ms1_isotopes = 2;
    /// Treat a non-CCS spectrum ion-mobility array as 1/K0, as in SearchOptions.
    bool ion_mobility_array_as_vssc = false;
  };

  /**
    @brief Stream spectra and keep the raw m/z, RT, 1/K0 and intensity points of selected precursors.

    MS1 points are kept around the precursor monoisotope and isotopes. MS2 points are kept around every
    library fragment, only from spectra whose isolation window contains the precursor m/z. A target is
    emitted to the sink as soon as the run passes its RT window, so memory is bounded by the targets that
    overlap the current RT. Points carry ion mobility -1 when the spectrum has no 1/K0.

    This is a read-only inspection side channel: it does not change any search result.
  */
  class PeakMapCollector final : public Interfaces::IMSDataConsumer
  {
  public:
    /// Receives one peak map and its MS level (1 precursor, 2 fragment).
    using Sink = std::function<void(const PeakMapExtractor::ExtractedPeakMap&, int)>;

    PeakMapCollector(std::vector<PeakMapTarget> targets, PeakMapOptions options, Sink sink);

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram&) override;
    void setExpectedSize(std::size_t, std::size_t) override;
    void setExperimentalSettings(const ExperimentalSettings&) override;

    /// Emit every target that is still open. Further spectra are ignored.
    void finish();

    /// Spectra whose RT was below an earlier spectrum. Targets already emitted may miss their points.
    std::size_t outOfOrderSpectra() const;
    /// Number of raw points passed to the sink.
    std::size_t pointsWritten() const;
    /// Number of peak maps passed to the sink.
    std::size_t peakMapsWritten() const;

  private:
    struct State
    {
      double rt_start = -std::numeric_limits<double>::infinity();
      double rt_end = std::numeric_limits<double>::infinity();
      std::vector<PeakMapExtractor::ExtractedPeakMap> ms1;
      std::vector<PeakMapExtractor::ExtractedPeakMap> ms2;
      bool emitted = false;
    };

    void extract_(const MSSpectrum& spectrum, double rt);
    void emit_(std::size_t target);

    std::vector<PeakMapTarget> targets_;
    PeakMapOptions options_;
    Sink sink_;
    std::vector<State> states_;
    std::vector<std::size_t> open_order_;
    std::size_t next_open_ = 0;
    std::vector<std::size_t> active_;
    double max_rt_ = -std::numeric_limits<double>::infinity();
    std::size_t out_of_order_ = 0;
    std::size_t points_written_ = 0;
    std::size_t peak_maps_written_ = 0;
    bool finished_ = false;
  };

  /// Forward every spectrum to two consumers in order (the search first).
  class TeeConsumer final : public Interfaces::IMSDataConsumer
  {
  public:
    TeeConsumer(Interfaces::IMSDataConsumer& first, Interfaces::IMSDataConsumer& second);

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram& chromatogram) override;
    void setExpectedSize(std::size_t spectra, std::size_t chromatograms) override;
    void setExperimentalSettings(const ExperimentalSettings& settings) override;

  private:
    Interfaces::IMSDataConsumer& first_;
    Interfaces::IMSDataConsumer& second_;
  };
}

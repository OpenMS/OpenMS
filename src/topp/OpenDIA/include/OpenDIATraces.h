// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /// A caller-labelled candidate window. Unknown targets must keep label -1.
  struct TraceWindow
  {
    std::string sample_id;
    std::string precursor_id;
    std::string group_id;
    int label = -1;
    double apex_rt = 0.0;
  };

  /**
    @brief Optional per-window extraction evidence for the experimental second-stage score.

    Collected only after TraceCollector::enableEvidence(). RT bins match the raw tensor.
    Mass errors are intensity-weighted signed ppm sums over points within the core window around the apex.
  */
  struct TraceEvidence
  {
    static constexpr std::size_t RT_BINS = 32;
    std::array<float, RT_BINS> ms1{};  ///< MS1 monoisotope intensity per RT bin.
    std::uint32_t ms1_mask = 0;        ///< RT bins that received an MS1 spectrum.
    std::uint32_t ms2_mask = 0;        ///< RT bins that received an MS2 spectrum isolating the precursor.
    float ms1_isotope = 0.0f;          ///< MS1 first-isotope intensity in the core window.
    double ms1_ppm_sum = 0.0;          ///< Sum of intensity x signed ppm, MS1 monoisotope, core window.
    double ms1_weight = 0.0;           ///< Sum of MS1 monoisotope intensity, core window.
    double ms2_ppm_sum = 0.0;          ///< Sum of intensity x signed ppm over fragments, core window.
    double ms2_weight = 0.0;           ///< Sum of fragment intensity, core window.
  };

  /**
    @brief Extract bounded raw RT/IM intensity tensors for a specified sample of candidates.

    This second streaming pass retains no run-wide spectra. Six library fragments
    are ordered by decreasing library intensity and binned into 32 RT by 8 IM cells.
    No peak-shape, coelution or similarity scores are supplied to the model.
  */
  class TraceCollector final : public Interfaces::IMSDataConsumer
  {
  public:
    static constexpr std::size_t FRAGMENTS = 6;
    static constexpr std::size_t RT_BINS = 32;
    static constexpr std::size_t IM_BINS = 8;
    static constexpr std::size_t SAMPLE_SIZE = FRAGMENTS * RT_BINS * IM_BINS;
    /// Default raw tensor budget (2 GiB).
    static constexpr std::size_t DEFAULT_MAX_TENSOR_BYTES = std::size_t{2} << 30;

    /** @brief Construct a bounded candidate tensor collector.
      @param[in] windows Explicit candidate windows, labels and related-peptide group IDs.
      @param[in] library Prepared transition library.
      @param[in] mz_tolerance Matching half width in Da.
      @param[in] rt_half_width RT half width in seconds.
      @param[in] im_half_width Mobility half width in 1/K0.
      @param[in] max_samples Maximum allowed sample count.
      @param[in] max_tensor_bytes Hard cap on the retained raw tensor size in bytes.
    */
    TraceCollector(std::vector<TraceWindow> windows, const OpenSwath::LightTargetedExperiment& library,
                   double mz_tolerance = 0.02, double rt_half_width = 15.0,
                   double im_half_width = 0.03, std::size_t max_samples = 16384,
                   std::size_t max_tensor_bytes = DEFAULT_MAX_TENSOR_BYTES);

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram&) override;
    void setExpectedSize(std::size_t, std::size_t) override;
    void setExperimentalSettings(const ExperimentalSettings&) override;

    /** @brief Write little-endian float32 NumPy tensors and matching sample metadata.
      @param[in] prefix Output prefix for .traces.npy, .samples.tsv and .geometry.json.
    */
    void store(const std::string& prefix) const;

    /// Raw tensor in sample, fragment, RT, mobility order, also consumed directly by the native CNN.
    const std::vector<float>& intensities() const { return intensities_; }
    /// Whether sample i contains compatible, library-centred mobility measurements.
    bool hasIM(std::size_t i) const { return samples_.at(i).has_im; }
    /// Candidate identity and explicit seed label, never network inputs.
    const TraceWindow& window(std::size_t i) const { return samples_.at(i).window; }
    /// Library decoy status, used only for validating seed labels.
    bool isDecoy(std::size_t i) const { return samples_.at(i).is_decoy; }
    /// Unnormalized library fragment intensities in tensor channel order.
    const std::array<double, FRAGMENTS>& libraryIntensities(std::size_t i) const { return samples_.at(i).library_intensity; }
    /// Fragment m/z in tensor channel order, with zero denoting padding.
    const std::array<double, FRAGMENTS>& fragmentMZ(std::size_t i) const { return samples_.at(i).fragment_mz; }
    /// Treat a spectrum ion-mobility array with a non-CCS unit as 1/K0 (VSSC).
    void setIonMobilityArrayAsVSSC(bool value) { ion_mobility_array_as_vssc_ = value; }

    /**
      @brief Also collect TraceEvidence (MS1 traces and m/z errors). Off by default; the tensor is unchanged.

      MS1 spectra must then be passed to consumeSpectrum(); without this call they are ignored.
      @param[in] ms1_ppm Half width of the MS1 m/z windows in ppm.
      @param[in] core_half_width_s Half width around the apex used for the mass errors and the isotope.
    */
    void enableEvidence(double ms1_ppm, double core_half_width_s);
    /// Whether enableEvidence() was called.
    bool evidenceEnabled() const { return !evidence_.empty(); }
    /// Evidence of sample i; requires evidenceEnabled().
    const TraceEvidence& evidence(std::size_t i) const { return evidence_.at(i); }
    /// RT half width of the tensor in seconds.
    double rtHalfWidth() const { return rt_half_width_; }

  private:
    struct Sample
    {
      TraceWindow window;
      double precursor_mz = 0.0;
      double precursor_im = 0.0;
      int precursor_charge = 0;
      bool is_decoy = false;
      bool has_im = false;
      bool saw_missing_im = false;
      std::array<double, FRAGMENTS> fragment_mz{};
      std::array<double, FRAGMENTS> library_intensity{};
    };
    std::vector<Sample> samples_;
    std::vector<std::size_t> by_rt_;
    std::vector<float> intensities_;
    double mz_tolerance_;
    double rt_half_width_;
    double im_half_width_;
    bool ion_mobility_array_as_vssc_ = false;
    std::vector<TraceEvidence> evidence_;
    double ms1_ppm_ = 0.0;
    double core_half_width_ = 0.0;

    void consumeMS1_(MSSpectrum& spectrum);
    std::size_t rtBin_(double rt, double apex_rt) const;
  };

  /** @brief Read explicitly selected windows from a TSV, rejecting missing grouping/labels.
    @param[in] path TSV containing sample_id, precursor_id, group_id, label, rt_apex.
    @param[in] max_samples Maximum rows to accept, not a silent truncation limit.
    @return Requested candidate windows in file order.
  */
  std::vector<TraceWindow> readTraceWindows(const std::string& path, std::size_t max_samples);
}

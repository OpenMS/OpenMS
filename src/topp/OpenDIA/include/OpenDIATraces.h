// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenDIALibrary.h>

#include <array>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
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
    /// With TraceCollector::enableByIonCounts(): |RT - apex| of the isolating MS2 spectrum the counts below come from.
    float byion_apex_dt = std::numeric_limits<float>::infinity();
    std::uint16_t byion_b = 0;         ///< Matched singly charged b ions in the isolating MS2 spectrum nearest the apex.
    std::uint16_t byion_y = 0;         ///< Matched singly charged y ions in the same spectrum.
  };

  /**
    @brief Optional per-window fragment spectrum evidence (TraceCollector::enablePeakGroupEvidence()).

    Sums over the isolating MS2 spectra within the core half width of the apex, points within the tolerance and the
    tensor's mobility window, per tensor channel: the monoisotope, its isotopes +1..+4 at the library fragment charge and
    the positions one 13C spacing before it at charges 1..4.
  */
  struct PeakGroupEvidence
  {
    static constexpr std::size_t FRAGMENTS = 6;
    static constexpr std::size_t ISOTOPES = 4;
    std::array<float, FRAGMENTS> mono{};    ///< Monoisotope intensity.
    std::array<float, FRAGMENTS> ppm_sum{}; ///< Sum of intensity x signed ppm of the monoisotope points.
    std::array<std::array<float, ISOTOPES>, FRAGMENTS> isotope{}; ///< Isotope i + 1 intensity.
    std::array<std::array<float, ISOTOPES>, FRAGMENTS> before{};  ///< Intensity at mono - 1.0033548 / (c + 1).
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
    /**
      @brief Streaming mode callback: receives windows that can no longer receive spectra (apex + RT half width before
      the current spectrum), while their tensors are still held. windowTensor(), evidence() and the per-window accessors
      are valid for these windows during the call; the tensors are released after it returns.
    */
    using CloseCallback = std::function<void(const std::vector<std::uint32_t>& closed, const TraceCollector& traces)>;

    /** @brief Construct a bounded candidate tensor collector.
      @param[in] windows Explicit candidate windows, labels and related-peptide group IDs.
      @param[in] library Prepared transition library.
      @param[in] mz_tolerance Matching half width in Da.
      @param[in] rt_half_width RT half width in seconds.
      @param[in] im_half_width Mobility half width in 1/K0.
      @param[in] max_samples Maximum allowed sample count.
      @param[in] max_tensor_bytes Hard cap on the retained raw tensor size in bytes (all windows; in streaming mode the
                 windows held at one time).
      @param[in] on_close If set, streaming mode: a window's tensor is held only from its first spectrum until it is
                 closed, then passed to @p on_close in batches and released, so memory follows the windows co-eluting
                 at one time instead of all windows. Spectra must arrive in RT order (an out-of-order spectrum that would
                 reach a closed window throws). Call finish() after the last spectrum. intensities() is empty.
    */
    TraceCollector(std::vector<TraceWindow> windows, const OpenDIA::Library& library,
                   double mz_tolerance = 0.02, double rt_half_width = 15.0,
                   double im_half_width = 0.03, std::size_t max_samples = 16384,
                   std::size_t max_tensor_bytes = DEFAULT_MAX_TENSOR_BYTES, CloseCallback on_close = {});

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram&) override;
    void setExpectedSize(std::size_t, std::size_t) override;
    void setExperimentalSettings(const ExperimentalSettings&) override;
    /// Streaming mode: close and pass on every remaining window. No effect otherwise.
    void finish();
    /// Streaming mode: largest number of tensors held at one time.
    std::size_t peakHeldTensors() const { return peak_held_; }

    /**
      @brief Stored raw tensor in sample, fragment, RT, mobility order, also consumed directly by the native CNN.

      Holds imBins() mobility bins per RT bin. When no candidate precursor has a library ion mobility, every point falls
      into the central of the 8 mobility bins, so only that bin is stored (imBins() == 1, [samples,6,32]). Use
      windowTensor() for the full [6,32,8] layout of one window.
    */
    const std::vector<float>& intensities() const { return intensities_; }
    /// Number of candidate windows.
    std::size_t size() const { return samples_.size(); }
    /// Stored mobility bins per RT bin: IM_BINS, or 1 when no candidate precursor has a library mobility (see intensities()).
    std::size_t imBins() const { return im_bins_; }
    /**
      @brief Window i in the full [6,32,8] layout.
      @param[in] i Window index.
      @param[out] buffer Scratch space, filled only when imBins() == 1.
      @return A pointer into the stored tensor, or into @p buffer.
    */
    const float* windowTensor(std::size_t i, std::array<float, SAMPLE_SIZE>& buffer) const;
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
    /// Library precursor m/z of sample i.
    double precursorMZ(std::size_t i) const { return samples_.at(i).precursor_mz; }
    /// Library precursor charge of sample i (0 if unknown).
    int precursorCharge(std::size_t i) const { return samples_.at(i).precursor_charge; }
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
    /**
      @brief Also count matched b and y ions (fragment charge 1) in the isolating MS2 spectrum nearest each window apex.

      Requires enableEvidence() first; only spectra within its core half width of the apex are considered. For every
      theoretical singly charged b and y ion of the library sequence (modifications included) between 150.06 and 2000 Th,
      the nearest spectrum peak within @p ppm counts as a match if its intensity exceeds @p intensity_min. Ion mobility is
      ignored (points are pooled). Targets and decoys are treated identically. Results are in TraceEvidence::byion_b/_y.

      @param[in] library The library passed to the constructor (sequences are read from its compounds).
      @param[in] ppm Matching half width in ppm.
      @param[in] intensity_min Minimum intensity of the matched peak.
      @throws Exception::InvalidParameter If evidence is not enabled, a tolerance is not positive or a sequence cannot be parsed.
    */
    void enableByIonCounts(const OpenDIA::Library& library, double ppm, double intensity_min);
    /// Evidence of sample i; requires evidenceEnabled().
    const TraceEvidence& evidence(std::size_t i) const { return evidence_.at(i); }
    /**
      @brief Also collect PeakGroupEvidence. Requires enableEvidence() first (its core half width is used).
      @param[in] ppm Matching half width in ppm around each target m/z.
      @throws Exception::InvalidParameter If evidence is not enabled or @p ppm is not positive.
    */
    void enablePeakGroupEvidence(double ppm);
    /// Peak-group evidence of sample i; requires enablePeakGroupEvidence().
    const PeakGroupEvidence& peakGroupEvidence(std::size_t i) const { return pg_evidence_.at(i); }
    /// Library fragment charge per tensor channel (0 if unknown or padding).
    const std::array<std::int8_t, FRAGMENTS>& fragmentCharges(std::size_t i) const { return samples_.at(i).fragment_charge; }
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
      std::array<std::int8_t, FRAGMENTS> fragment_charge{};
    };
    std::vector<Sample> samples_;
    std::vector<std::size_t> by_rt_;
    /// Samples by precursor m/z bucket (MZ_BUCKET_TH wide, from mz_bucket_origin_), each bucket in by_rt_ order.
    static constexpr double MZ_BUCKET_TH = 1.0;
    std::vector<std::vector<std::uint32_t>> by_mz_rt_;
    double mz_bucket_origin_ = 0.0;
    std::vector<float> intensities_;
    std::size_t im_bins_ = IM_BINS;
    double mz_tolerance_;
    double rt_half_width_;
    double im_half_width_;
    bool ion_mobility_array_as_vssc_ = false;
    std::vector<TraceEvidence> evidence_;
    double ms1_ppm_ = 0.0;
    double core_half_width_ = 0.0;
    std::vector<float> byion_mz_;                           ///< Pooled ion m/z: per sample, b ions then y ions.
    std::vector<std::array<std::uint32_t, 3>> byion_range_; ///< Per sample: first b, first y, end in byion_mz_.
    double byion_ppm_ = 0.0;
    double byion_intensity_min_ = 0.0;
    std::vector<PeakGroupEvidence> pg_evidence_;
    double pg_ppm_ = 0.0;
    // Streaming mode (see the constructor).
    static constexpr std::uint32_t NO_SLOT = std::numeric_limits<std::uint32_t>::max();
    static constexpr std::size_t SLOT_CHUNK = 16384;
    static constexpr std::size_t CLOSE_BATCH = 4096;
    CloseCallback on_close_;
    std::size_t stride_ = 0;                         ///< Stored floats per window.
    std::size_t max_held_ = 0;                       ///< Tensor budget in windows (streaming mode).
    std::vector<std::unique_ptr<float[]>> chunks_;   ///< Slot storage, SLOT_CHUNK slots per chunk.
    std::vector<std::uint32_t> slot_of_;
    std::vector<std::uint32_t> free_slots_;
    std::vector<std::uint8_t> closed_;
    std::vector<std::uint32_t> pending_;
    std::size_t close_cursor_ = 0;
    std::size_t held_ = 0;
    std::size_t peak_held_ = 0;

    float* tensor_(std::size_t index);
    const float* storedTensor_(std::size_t index) const;
    void close_(double rt);
    void flush_();

    void consumeMS1_(MSSpectrum& spectrum);
    void countByIons_(const MSSpectrum& spectrum, std::size_t index, double dt);
    std::size_t rtBin_(double rt, double apex_rt) const;
  };

  /** @brief Read explicitly selected windows from a TSV, rejecting missing grouping/labels.
    @param[in] path TSV containing sample_id, precursor_id, group_id, label, rt_apex.
    @param[in] max_samples Maximum rows to accept, not a silent truncation limit.
    @return Requested candidate windows in file order.
  */
  std::vector<TraceWindow> readTraceWindows(const std::string& path, std::size_t max_samples);
}
